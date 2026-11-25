# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "pyyaml"]
# ///
"""Test resources: OP3.

The Open Problems Perturbation Prediction (OP3) dataset: PBMCs of several donors treated with
small molecules, profiled with scRNA-seq (GEO GSE279945). It is the largest single-cell small
molecule perturbation dataset in primary tissue with multiple donor replicates.

Run from the root of the repository:

    uv run scripts/create_test_resources/openproblems_op3.py
"""

import anndata as ad
import lamindb as ln

import helpers

helpers.check_instance()
ln.track("gxI8Ed6sqoic", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/openproblems_op3.py")

# load the dataset from GEO and subsample it to 600 cells and 1500 genes; the loader reads the
# 14 GB h5ad into memory
output_dir = helpers.run_workflow("scrnaseq/load_openproblems_op3", {
  "id": "openproblems_op3",
  "input": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad",
  "dataset_name": "OP3: single-cell multimodal dataset in PBMCs for perturbation prediction benchmarking",
  "dataset_summary": "The Open Problems Perturbation Prediction (OP3) dataset with small molecule perturbations in PBMCs",
  "dataset_description": (
    "The OP3 dataset is to-date the largest single-cell small molecule perturbation dataset in "
    "primary tissue with multiple donor replicates."
  ),
  "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945",
  "dataset_reference": "GSE279945",
  "do_subsample": True,
  "n_obs": 600,
  "n_vars": 1500,
  "seed": 123,
  "output_raw": "$id/raw.anndata.zarr",
  "output_state": "$id/load_state.yaml",
}, highmem="100.GB")

output_dir = helpers.run_workflow("scrnaseq/process_scrnaseq", {
  "id": "openproblems_op3",
  "input": str(output_dir / "openproblems_op3" / "raw.anndata.zarr"),
  "output_dataset": "$id/dataset.anndata.zarr",
  "output_state": "$id/process_state.yaml",
}, output_dir=output_dir)

print(ad.read_zarr(output_dir / "openproblems_op3" / "dataset.anndata.zarr"))

# store the raw and the processed dataset as artifacts, in a collection
helpers.save_collection(
  output_dir / "openproblems_op3",
  ["raw.anndata.zarr", "dataset.anndata.zarr"],
  key="resources_test/datasets/openproblems_op3",
  modality="scrnaseq",
  description="OP3 (PBMCs, small molecule perturbations), subsampled to 600 cells and 1500 genes",
)

ln.finish()
