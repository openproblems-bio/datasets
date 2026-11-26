"""Tests for annotate_artifacts.py that don't need a LaminDB instance.

    uv run --with pytest --with "lamindb[aws,bionty]==2.12.2" --with anndata==0.13.2 \
      --with mudata==0.4.1 --with pyyaml pytest scripts/lamin/test_annotate_artifacts.py
"""

import sys
from pathlib import Path

import bionty as bt
import pytest

sys.path.insert(0, str(Path(__file__).parent))
import annotate_artifacts  # noqa: E402


def test_describe_dataframe_with_nullable_string_index(tmp_path):
  # with pandas 3, anndata writes string indexes as nullable string arrays (a group), as
  # the OP3 loader does
  import anndata as ad
  import numpy as np
  import pandas as pd
  import zarr

  obs = pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}, index=pd.Index(["c1", "c2", None], dtype="string"))
  with ad.settings.override(allow_write_nullable_strings=True):
    ad.AnnData(X=np.zeros((3, 2)), obs=obs).write_zarr(tmp_path / "test.zarr")
  root = zarr.open_group(tmp_path / "test.zarr", mode="r")
  assert isinstance(root["obs"][root["obs"].attrs["_index"]], zarr.Group)

  assert annotate_artifacts._describe_element(root["obs"]) == {"type": "dataframe", "shape": [3, 2]}


@pytest.mark.parametrize("current", [
  None,
  # set by the output rules in src/config/lamin_denbi.config
  "Raw scRNA-seq dataset",
  "Common multimodal dataset",
  # nf-lamin's default
  "Output artifact for run 12",
])
def test_generic_description_becomes_summary(current):
  assert annotate_artifacts._description(current, "A summary") == "A summary"


def test_specific_description_is_kept():
  assert annotate_artifacts._description("Immune Cell Atlas, three donors", "A summary") == "Immune Cell Atlas, three donors"
  assert annotate_artifacts._description(None, None) is None


def test_slot_schema_name():
  items = [{"name": "normalized"}, {"name": "counts"}]
  assert annotate_artifacts.slot_schema_name("layers", items) == "layers: counts, normalized"
  items = [{"name": f"column_{i}"} for i in range(5)]
  assert annotate_artifacts.slot_schema_name("obs", items) == "obs: 5 features"


def _content(slots, keys=None):
  # what lamindb hashes a composite schema by: its slot keys and their features, not the names
  return frozenset(
    (key, slot, frozenset((i["name"], i["type"], i.get("multiple", False), i.get("required", True)) for i in items))
    for key, (slot, items) in slots.items()
    if keys is None or key in keys
  )


def test_formats_have_distinct_composite_schemas():
  # a composite schema is named after its format, so two formats with the same slots would
  # share one schema, renamed to whichever format was annotated last
  repo = Path(__file__).parents[2]
  full, curated = {}, {}
  for path in annotate_artifacts.SPECS.values():
    fmt = annotate_artifacts.load_spec(repo / path)["info"]["format"]
    slots = annotate_artifacts.slot_items(fmt)
    full.setdefault(_content(slots), []).append(path)
    for modalities in [["rna", "prot"], ["rna", "atac"]] if fmt["type"] == "mudata_zarr" else [None]:
      keys = annotate_artifacts.curated_slot_names(modalities)
      curated.setdefault(_content(slots, keys), []).append((path, modalities))
  assert [paths for paths in full.values() if len(paths) > 1] == []
  assert [paths for paths in curated.values() if len(paths) > 1] == []


def test_placeholder_registries_cover_cellxgene_unknowns():
  # CELLxGENE schema 5 allows "unknown" in these columns (and "na" for ethnicity), so they
  # need a placeholder term to validate; the Census release of 2025-11-08 has unknown cell types
  assert annotate_artifacts.PLACEHOLDER_REGISTRIES == {
    "cell_type_ontology_term_id": bt.CellType,
    "development_stage_ontology_term_id": bt.DevelopmentalStage,
    "self_reported_ethnicity_ontology_term_id": bt.Ethnicity,
    "sex_ontology_term_id": bt.Phenotype,
  }
