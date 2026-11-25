# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws]==2.12.2", "pyyaml"]
# ///
"""Process the latest raw datasets of a modality in LaminDB on de.NBI.

Finds the latest version of every `<modality>/<id>/raw.{anndata,mudata}.zarr` artifact in the
instance and launches `process_<modality>` on them, reading each raw dataset by its uid so
LaminDB links the processed dataset to it. Datasets marked `large` in the load scripts get
more memory and time.

Run from the root of the repository, after the load runs have finished:

    uv run scripts/create_resources/process_datasets.py scrnaseq [--ids ...] [--dry-run]

and annotate the results once the runs have finished:

    uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --key-prefix scrnaseq/
"""

import re

import cellxgene_census
import helpers
import openproblems_neurips2021_bmmc
import openproblems_neurips2022_pbmc
import openproblems_op3

LARGE_IDS = {
  dataset["id"]
  for module in [cellxgene_census, openproblems_op3, openproblems_neurips2021_bmmc, openproblems_neurips2022_pbmc]
  for dataset in module.DATASETS
  if dataset.get("large")
}

RAW_KEY = re.compile(r"(?P<modality>scrnaseq|multimodal)/(?P<id>.+)/raw\.(?P<format>anndata|mudata)\.zarr")


def raw_key_to_param(key, uid, instance):
  """The process parameters of the raw dataset artifact with this key and uid."""
  match = RAW_KEY.fullmatch(key)
  if match is None:
    raise ValueError(f"Not the key of a raw dataset: '{key}'")
  return {
    "id": match["id"],
    "input": f"lamin://{instance}/artifact/{uid}",
    "output_dataset": f"$id/dataset.{match['format']}.zarr",
  }


# Viash reads the state file name from the top-level parameters only
SHARED = {"output_state": "$id/process_state.yaml"}


def latest_raw_datasets(instance, modality):
  import lamindb as ln

  datasets = []
  for artifact in ln.Artifact.connect(instance).filter(key__startswith=f"{modality}/", is_latest=True):
    try:
      datasets.append(raw_key_to_param(artifact.key, artifact.uid, instance))
    except ValueError:
      continue
  return sorted(datasets, key=lambda dataset: dataset["id"])


def main():
  parser = helpers.make_parser(__doc__)
  parser.add_argument("modality", choices=helpers.PREFIXES)
  args = parser.parse_args()

  instance = helpers.resolve_instance(args.instance)
  datasets = [
    {**dataset, "large": dataset["id"] in LARGE_IDS}
    for dataset in latest_raw_datasets(instance, args.modality)
  ]
  if not datasets:
    raise SystemExit(f"No raw {args.modality} datasets in {instance}")
  helpers.launch_datasets(f"{args.modality}/process_{args.modality}", datasets, SHARED, args.modality, args)
  print(f"\nOnce the runs have finished:\n  uv run scripts/lamin/annotate_artifacts.py --instance {instance} --key-prefix {args.modality}/")


if __name__ == "__main__":
  main()
