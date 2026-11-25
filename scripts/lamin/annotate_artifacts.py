# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "lamindb[aws,bionty]==2.12.2",
#   "anndata==0.13.2",
#   "mudata==0.4.1",
#   "pyyaml",
# ]
# ///
"""Annotate LaminDB artifacts with everything the OpenProblems website shows about a dataset.

Run it after anything that published datasets to LaminDB, be it a Nextflow run or a script:

    uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --run <run uid>
    uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --key-prefix scrnaseq/
    uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --artifact <uid> --spec src/api/file_raw.yaml

Validation stays in OpenProblems (`validate_dataset`); this records the result in LaminDB, so
the hub page and the website can be rendered from the instance alone:

* a schema generated from the `src/api` file format, with every slot (`obs`, `var`, `uns`,
  `layers`, `obsm`, ...), so the structure of the file is documented next to the data
* ontology terms of the `*_ontology_term_id` columns, linked through bionty
* the `dataset_*` metadata and processing parameters as feature values
* the actual structure of the store (slot, name, dtype, shape) and cell counts per category
* references from the bibliography, with their DOI, and the modality

For an AnnData store, only metadata, `obs`, `var` and `uns` are read, straight from where the
artifact is stored, so this works for datasets of any size. lamindb loads a MuData store into
memory to validate it. Annotating is idempotent: artifacts that already have the schema of
their file format are skipped.

Each file format has its own reader for the structure and summary of a store: `anndata_zarr`
and `mudata_zarr` exist, SpatialData, data frames and JSON/YAML get one later.
"""

import argparse
import re
from pathlib import Path

import anndata as ad
import bionty as bt
import lamindb as ln
import yaml
import zarr

REPO_ROOT = Path(__file__).resolve().parents[2]
BIB_PATH = REPO_ROOT / "common/library.bib"
PROJECT = "openproblems-bio/datasets"

# file format of each file, by the last part of its key
SPECS = {
  "raw.anndata.zarr": "src/api/file_raw.yaml",
  "normalized.anndata.zarr": "src/api/file_normalized.yaml",
  "hvg.anndata.zarr": "src/api/file_hvg.yaml",
  "pca.anndata.zarr": "src/api/file_pca.yaml",
  "dataset.anndata.zarr": "src/api/file_common_dataset.yaml",
  "raw.mudata.zarr": "src/api/file_multimodal_raw.yaml",
  "dataset.mudata.zarr": "src/api/file_multimodal_dataset.yaml",
}

# modalities, by the first part of a key
MODALITIES = ["scrnaseq", "multimodal", "spatial"]

# spec types to lamindb dtypes
DTYPES = {
  "string": str,
  "categorical": str,
  "integer": int,
  "long": int,
  "double": float,
  "float": float,
  "boolean": bool,
  "object": dict,
}

# obs columns validated against an ontology, as in the CELLxGENE schema
ONTOLOGY_COLUMNS = {
  "assay_ontology_term_id": "cat[bionty.ExperimentalFactor.ontology_id]",
  "cell_type_ontology_term_id": "cat[bionty.CellType.ontology_id]",
  "development_stage_ontology_term_id": "cat[bionty.DevelopmentalStage.ontology_id]",
  # CELLxGENE uses PATO:0000461 ("normal") next to MONDO diseases
  "disease_ontology_term_id": "cat[bionty.Disease.ontology_id|bionty.Phenotype.ontology_id]",
  "organism_ontology_term_id": "cat[bionty.Organism.ontology_id]",
  "self_reported_ethnicity_ontology_term_id": "cat[bionty.Ethnicity.ontology_id]",
  "sex_ontology_term_id": "cat[bionty.Phenotype.ontology_id]",
  "tissue_ontology_term_id": "cat[bionty.Tissue.ontology_id]",
  "tissue_general_ontology_term_id": "cat[bionty.Tissue.ontology_id]",
}

# registries in which to create the placeholder terms CELLxGENE uses ("unknown", "na")
PLACEHOLDER_REGISTRIES = {
  "cell_type_ontology_term_id": bt.CellType,
  "development_stage_ontology_term_id": bt.DevelopmentalStage,
  "self_reported_ethnicity_ontology_term_id": bt.Ethnicity,
  "sex_ontology_term_id": bt.Phenotype,
}

# obs columns of which the number of cells per category is stored
COUNT_COLUMNS = [
  "assay", "cell_type", "development_stage", "disease", "donor_id", "organism",
  "self_reported_ethnicity", "sex", "tissue", "tissue_general", "batch",
]

CURATED_SLOTS = ["obs", "var", "uns"]
ANNDATA_SLOTS = ["X", "obs", "var", "layers", "obsm", "obsp", "varm", "varp", "uns"]
# slots of the MuData object itself; its var is the concatenation of the modalities' var
MUDATA_SLOTS = ["obs", "obsm", "obsp", "uns"]
ROOT_TYPE = "OpenProblems"


def load_spec(path):
  """Read a `src/api` file spec, resolving `__merge__` like viash does."""
  path = Path(path)
  spec = yaml.safe_load(path.read_text())
  parent = spec.pop("__merge__", None)
  if parent is None:
    return spec
  return _merge(load_spec((path.parent / parent).resolve()), spec)


def _merge(a, b):
  if isinstance(a, dict) and isinstance(b, dict):
    return {**a, **{k: _merge(a[k], v) if k in a else v for k, v in b.items()}}
  if isinstance(a, list) and isinstance(b, list):
    return a + b
  return b


# --- registries ----------------------------------------------------------------------------

def _get_or_create(registry, **fields):
  record = registry.filter(**fields).one_or_none()
  return record if record is not None else registry(**fields).save()


def _feature_type(name):
  """Features are grouped per slot, since e.g. `dataset_id` means something else in obs and uns."""
  root = _get_or_create(ln.Feature, name=ROOT_TYPE, is_type=True)
  return _get_or_create(ln.Feature, name=name, is_type=True, type=root)


def _feature(name, dtype, type_name, description=None):
  ftype = _feature_type(type_name)
  feature = ln.Feature.filter(name=name, type=ftype).one_or_none()
  if feature is None:
    feature = ln.Feature(name=name, dtype=dtype, type=ftype, description=description).save()
  return feature


def _spec_feature(item, slot):
  name = item["name"]
  if slot == "obs" and name in ONTOLOGY_COLUMNS:
    dtype = ONTOLOGY_COLUMNS[name]
  elif slot == "uns" and item.get("multiple"):
    dtype = list[DTYPES[item["type"]]]
  else:
    dtype = DTYPES[item["type"]]
  feature = _feature(name, dtype, slot, (item.get("description") or "").strip() or None)
  return feature if item.get("required", True) else feature.with_config(optional=True)


def _slot_schemas(name, fmt, slot_names, prefix=""):
  """One schema per slot of a format, named `<prefix><slot>` as lamindb's curators expect."""
  return {
    f"{prefix}{slot}": ln.Schema(
      name=f"{name}: {prefix}{slot}",
      features=[_spec_feature(item, slot) for item in items],
      otype="DataFrame" if slot in ("obs", "var") else None,
    ).save()
    for slot, items in fmt.items()
    if slot in slot_names and isinstance(items, list) and items
  }


def spec_schemas(spec_path, modalities=None):
  """Turn a file spec into a LaminDB schema with every slot, and one with the slots lamindb can curate.

  For a MuData format, the slots of a modality are named `<modality>:<slot>`. Only the
  modalities in `modalities` are curated, since a dataset has a subset of the optional ones.
  """
  spec = load_spec(spec_path)
  fmt = spec["info"]["format"]
  name = spec.get("label") or Path(spec_path).stem
  description = (spec.get("summary") or "").strip() or None
  if fmt["type"] == "mudata_zarr":
    otype = "MuData"
    slots = _slot_schemas(name, fmt, MUDATA_SLOTS)
    for modality, mod_fmt in (fmt.get("mod") or {}).items():
      slots.update(_slot_schemas(name, mod_fmt, ANNDATA_SLOTS, prefix=f"{modality}:"))
    curated_names = CURATED_SLOTS + [f"{m}:{slot}" for m in (modalities or []) for slot in CURATED_SLOTS]
  else:
    otype = "AnnData"
    slots = _slot_schemas(name, fmt, ANNDATA_SLOTS)
    curated_names = CURATED_SLOTS
  full = ln.Schema(name=name, otype=otype, slots=slots, description=description).save()
  curated_slots = {slot: schema for slot, schema in slots.items() if slot in curated_names}
  suffix = f": {', '.join(modalities)}" if modalities else ""
  curated = ln.Schema(name=f"{name} (curated slots{suffix})", otype=otype, slots=curated_slots).save()
  return spec, full, curated


# --- reading a store -----------------------------------------------------------------------

def _open_store(artifact):
  """Open the zarr store of an artifact where it is stored, without downloading it."""
  path = artifact.path
  if getattr(path, "protocol", "") in ("", "file", "local"):
    return zarr.open_group(str(path), mode="r")
  return zarr.open_group(zarr.storage.FsspecStore.from_upath(path, read_only=True), mode="r")


def _describe_element(node):
  """Type, dtype and shape of an element in an anndata zarr store, from its metadata only."""
  if isinstance(node, zarr.Array):
    return {"type": node.attrs.get("encoding-type", "array"), "dtype": str(node.dtype), "shape": list(node.shape)}
  encoding = node.attrs.get("encoding-type")
  if encoding in ("csr_matrix", "csc_matrix"):
    return {"type": encoding, "dtype": str(node["data"].dtype), "shape": list(node.attrs["shape"])}
  if encoding == "categorical":
    return {
      "type": "categorical",
      "dtype": str(node["categories"].dtype),
      "shape": list(node["codes"].shape),
      "n_categories": node["categories"].shape[0],
    }
  if encoding and encoding.startswith("nullable-"):
    return {"type": encoding, "dtype": str(node["values"].dtype), "shape": list(node["values"].shape)}
  if encoding == "dataframe":
    # the index is an array, or a group for nullable strings and categoricals
    n_rows = _describe_element(node[node.attrs["_index"]])["shape"][0]
    return {"type": "dataframe", "shape": [n_rows, len(node.attrs["column-order"])]}
  return {"type": encoding or "dict", "keys": sorted(node.keys())}


def _anndata_structure(root, spec, slots=ANNDATA_SLOTS):
  """Per slot, the name, type, dtype and shape of each element, with the description from the spec."""
  fmt = spec["info"]["format"]
  structure = {}
  for slot in slots:
    if slot not in root:
      continue
    descriptions = {item["name"]: (item.get("description") or "").strip() for item in fmt.get(slot, []) or [] if isinstance(item, dict)}
    if slot == "X":
      structure["X"] = [{"name": "X", **_describe_element(root["X"])}]
      continue
    group = root[slot]
    names = group.attrs["column-order"] if slot in ("obs", "var") else sorted(group.keys())
    structure[slot] = [
      {"name": name, **_describe_element(group[name]), "description": descriptions.get(name)}
      for name in names
    ]
  return structure


def _anndata_summary(root, structure):
  """Dataset-level numbers and counts that are shown next to the data."""
  obs = ad.io.read_elem(root["obs"])
  var = ad.io.read_elem(root["var"])
  uns = ad.io.read_elem(root["uns"]) if "uns" in root else {}
  summary = {"n_vars": int(var.shape[0])}
  for column, feature in [("donor_id", "n_donors"), ("batch", "n_batches")]:
    if column in obs:
      summary[feature] = int(obs[column].nunique())
  if "hvg" in var:
    summary["n_hvg"] = int(var["hvg"].sum())
  obsm = {item["name"]: item for item in structure.get("obsm", [])}
  if "X_pca" in obsm:
    summary["n_pca_components"] = int(obsm["X_pca"]["shape"][1])
  knn_params = uns.get("knn", {}).get("params", {})
  if "n_neighbors" in knn_params:
    summary["knn_n_neighbors"] = int(knn_params["n_neighbors"])
  summary["obs_counts"] = {
    column: {str(k): int(v) for k, v in obs[column].value_counts().items()}
    for column in COUNT_COLUMNS
    if column in obs
  }
  return obs, summary, uns


def _mudata_structure(root, spec):
  """The structure of the MuData object itself, and of each modality."""
  mod_fmts = spec["info"]["format"].get("mod") or {}
  structure = _anndata_structure(root, spec, slots=MUDATA_SLOTS)
  structure["mod"] = {
    modality: _anndata_structure(root["mod"][modality], {"info": {"format": mod_fmts.get(modality) or {}}})
    for modality in root["mod"].attrs.get("mod-order", sorted(root["mod"].keys()))
  }
  return structure


def _mudata_summary(root, structure):
  """Dataset-level numbers from the global obs, and per modality its size and processing."""
  obs = ad.io.read_elem(root["obs"])
  uns = ad.io.read_elem(root["uns"]) if "uns" in root else {}
  modalities = {}
  for modality, mod_structure in structure["mod"].items():
    var = ad.io.read_elem(root["mod"][modality]["var"])
    mod_uns = ad.io.read_elem(root["mod"][modality]["uns"]) if "uns" in root["mod"][modality] else {}
    info = {"n_vars": int(var.shape[0])}
    if "hvg" in var:
      info["n_hvg"] = int(var["hvg"].sum())
    obsm = {item["name"]: item for item in mod_structure.get("obsm", [])}
    if "X_svd" in obsm:
      info["n_svd_components"] = int(obsm["X_svd"]["shape"][1])
    if "normalization_id" in mod_uns:
      info["normalization_id"] = str(mod_uns["normalization_id"])
    modalities[modality] = info
  summary = {"n_vars": sum(info["n_vars"] for info in modalities.values()), "modalities": modalities}
  for column, feature in [("donor_id", "n_donors"), ("batch", "n_batches")]:
    if column in obs:
      summary[feature] = int(obs[column].nunique())
  summary["obs_counts"] = {
    column: {str(k): int(v) for k, v in obs[column].value_counts().items()}
    for column in COUNT_COLUMNS
    if column in obs
  }
  return obs, summary, uns


def _anndata_modalities(root):
  return None


def _mudata_modalities(root):
  return list(root["mod"].attrs.get("mod-order", sorted(root["mod"].keys())))


FORMATS = {
  "anndata_zarr": {
    "structure": _anndata_structure,
    "summary": _anndata_summary,
    "modalities": _anndata_modalities,
    "curator": ln.curators.AnnDataCurator,
  },
  # lamindb loads a MuData artifact into memory to curate it, rather than streaming it
  "mudata_zarr": {
    "structure": _mudata_structure,
    "summary": _mudata_summary,
    "modalities": _mudata_modalities,
    "curator": ln.curators.MuDataCurator,
  },
}


# --- annotating ----------------------------------------------------------------------------

def _register_organisms(ontology_ids):
  """Register the species itself for each organism term.

  Ensembl lists mouse strains under the same NCBI Taxonomy term as mouse, so bionty can't
  pick a record, and lamindb then loops forever saving the duplicates.
  """
  for ontology_id in ontology_ids:
    if bt.Organism.filter(ontology_id=ontology_id).exists():
      continue
    records = bt.Organism.from_source(ontology_id=ontology_id)
    if isinstance(records, list):
      records = min(records, key=lambda record: len(record.name))
    records.save()


def _register_missing_terms(registry, column, values, obs):
  """Register terms that are not in the default ontology of a registry.

  Mouse development stages are in MmusDv rather than the default HsapDv, so other sources
  are tried first. What's left are CELLxGENE placeholders ("unknown", "na"), which get the
  label of the matching obs column as name.
  """
  label_column = column.removesuffix("_ontology_term_id")
  labels = dict(zip(obs[column].astype(str), obs[label_column].astype(str))) if label_column in obs else {}
  sources = bt.Source.filter(entity=f"bionty.{registry.__name__}").order_by("-currently_used")
  for value in values:
    record = None
    for source in sources:
      try:
        record = registry.from_source(ontology_id=value, source=source, mute=True)
      except Exception:
        record = None
      if record is not None and not isinstance(record, list):
        break
      record = None
    if record is None:
      record = registry(name=labels.get(value, value), ontology_id=value, description="CELLxGENE placeholder term")
    record.save()


def _curate(artifact, curated_schema, curator_class, obs):
  """Validate obs/var/uns with lamindb and link the ontology terms.

  The curator opens an artifact in remote storage instead of loading it.
  """
  if "organism_ontology_term_id" in obs:
    _register_organisms(obs["organism_ontology_term_id"].unique())
  curator = curator_class(artifact, curated_schema)
  try:
    curator.validate()
  except ln.errors.ValidationError:
    for column, values in curator.slots["obs"].cat.non_validated.items():
      registry = PLACEHOLDER_REGISTRIES.get(column)
      if registry is None:
        raise
      _register_missing_terms(registry, column, values, obs)
    curator.validate()
  curator.save_artifact()


def _bib_entry(bib_text, key):
  match = re.search(r"@\w+\s*\{\s*" + re.escape(key) + r"\s*,(.*?)\n\}", bib_text, flags=re.S)
  if match is None:
    return None
  fields = {}
  for field, value in re.findall(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", match.group(1)):
    fields[field.lower()] = re.sub(r"[{}]", "", " ".join(value.split()))
  return fields


def _link_references(artifact, keys, bib_path):
  bib_text = Path(bib_path).read_text()
  for key in keys:
    entry = _bib_entry(bib_text, key)
    if entry is None:
      print(f"  reference '{key}' not found in {bib_path}")
      continue
    doi = entry.get("doi")
    reference = ln.Reference.filter(doi=doi).one_or_none() if doi else None
    if reference is None:
      reference = ln.Reference.filter(abbr=key).one_or_none()
    if reference is None:
      reference = ln.Reference(
        name=entry.get("title", key),
        abbr=key,
        doi=doi,
        url=entry.get("url") or (f"https://doi.org/{doi}" if doi else None),
      ).save()
    artifact.references.add(reference)


def annotate_artifact(artifact, spec_path, modality=None, bib_path=BIB_PATH, force=False):
  """Record what OpenProblems knows about a dataset artifact in LaminDB.

  Returns whether the artifact was annotated; it is skipped when it already has the schema
  of its file format, unless `force` is set.
  """
  spec, full_schema, _ = spec_schemas(spec_path)
  fmt = FORMATS[spec["info"]["format"]["type"]]
  if artifact.schema_id == full_schema.id and not force:
    return False

  root = _open_store(artifact)
  _, _, curated_schema = spec_schemas(spec_path, modalities=fmt["modalities"](root))
  structure = fmt["structure"](root, spec)
  obs, summary, uns = fmt["summary"](root, structure)

  # validate the slots lamindb can curate, and link their ontology terms
  _curate(artifact, curated_schema, fmt["curator"], obs)

  # the metadata in uns, the summary and the structure as feature values
  uns_features = {feature.name: feature for feature in full_schema.slots["uns"].members} if "uns" in full_schema.slots else {}
  values = {}
  for name, value in uns.items():
    if name not in uns_features:
      continue
    if hasattr(value, "tolist"):
      value = value.tolist()
    if isinstance(value, (str, int, float, bool, list)):
      values[uns_features[name]] = value
  for name, value in {**summary, "structure": structure}.items():
    values[_feature(name, dict if isinstance(value, dict) else int, "dataset")] = value
  if force:
    _remove_values(artifact, values)
  artifact.features.add_values(values)

  # modality, for grouping datasets on the website
  if modality is not None:
    modality_type = _get_or_create(ln.ULabel, name="Modality", is_type=True)
    artifact.ulabels.add(_get_or_create(ln.ULabel, name=modality, type=modality_type))

  # references, with their DOI
  references = uns.get("dataset_reference")
  if references is not None:
    _link_references(artifact, [references] if isinstance(references, str) else list(references), bib_path)

  # document every slot of the file format, not only the curated ones; this comes last, as
  # an artifact with the schema of its file format counts as annotated
  artifact.schema = full_schema
  artifact.n_observations = int(obs.shape[0])
  if artifact.description is None and uns.get("dataset_summary"):
    artifact.description = uns["dataset_summary"]
  artifact.save()
  return True


def _remove_values(artifact, values):
  """Remove values of an earlier annotation, so they are replaced rather than added to."""
  existing = artifact.features.get_values()
  for feature in values:
    if feature.name not in existing:
      continue
    try:
      artifact.features.remove_values(feature)
    except ValueError:
      # lamindb doesn't remove values of features in the artifact's schema
      pass


# --- command line --------------------------------------------------------------------------

def spec_for(artifact):
  """The file format of an artifact, from the last part of its key."""
  spec = SPECS.get(artifact.key.rsplit("/", 1)[-1]) if artifact.key else None
  return REPO_ROOT / spec if spec else None


def modality_for(artifact):
  """The modality of an artifact, from the first part of its key."""
  first = artifact.key.split("/", 1)[0] if artifact.key else None
  return first if first in MODALITIES else None


def select_artifacts(run=None, key_prefix=None, uids=None):
  if run is not None:
    return list(ln.Artifact.filter(run__uid=run))
  if key_prefix is not None:
    return list(ln.Artifact.filter(key__startswith=key_prefix, is_latest=True))
  return [ln.Artifact.get(uid) for uid in uids]


def main():
  parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
  parser.add_argument("--instance", required=True, help="LaminDB instance, e.g. openproblems/data-dev")
  which = parser.add_mutually_exclusive_group(required=True)
  which.add_argument("--run", help="annotate the artifacts created by this run (uid)")
  which.add_argument("--key-prefix", help="annotate the latest artifacts with keys starting with this prefix")
  which.add_argument("--artifact", nargs="+", help="annotate these artifacts (uids)")
  parser.add_argument("--spec", help="file format to use instead of the one derived from the key")
  parser.add_argument("--modality", choices=MODALITIES, help="modality to use instead of the one derived from the key")
  parser.add_argument("--force", action="store_true", help="annotate artifacts that are already annotated")
  parser.add_argument("--dry-run", action="store_true", help="only list what would be annotated")
  args = parser.parse_args()

  # connect explicitly: an instance without local settings silently falls back to the default
  ln.connect(args.instance)
  if ln.setup.settings.instance.slug != args.instance:
    raise RuntimeError(f"Connected to '{ln.setup.settings.instance.slug}' instead of '{args.instance}'")

  artifacts = select_artifacts(run=args.run, key_prefix=args.key_prefix, uids=args.artifact)
  todo = []
  for artifact in artifacts:
    spec = Path(args.spec).resolve() if args.spec else spec_for(artifact)
    modality = args.modality or modality_for(artifact)
    if spec is None:
      print(f"skip {artifact.uid} {artifact.key}: no file format for this key")
      continue
    todo.append((artifact, spec, modality))
    print(f"{artifact.uid} {artifact.key}: {spec.relative_to(REPO_ROOT)}, modality {modality}")
  if args.dry_run or not todo:
    return

  ln.track("rduBFeXpb7eO", project=PROJECT)
  for artifact, spec, modality in todo:
    annotated = annotate_artifact(artifact, spec, modality=modality, force=args.force)
    print(f"{artifact.uid}: {'annotated' if annotated else 'already annotated'}")
  ln.finish()


if __name__ == "__main__":
  main()
