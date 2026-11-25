# /// script
# requires-python = ">=3.12"
# dependencies = ["pyyaml"]
# ///
"""Load the OP3 dataset on de.NBI and publish it to LaminDB.

The loader reads the 14 GB h5ad from GEO into memory, so it runs with more memory and time.

Run from the root of the repository:

    uv run scripts/create_resources/openproblems_op3.py [--dry-run]
"""

import helpers

DATASETS = [
  {
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
    "large": True,
  },
]

SHARED = {
  "output_raw": "$id/raw.anndata.zarr",
  "output_state": "$id/load_state.yaml",
}

if __name__ == "__main__":
  args = helpers.parse_args(__doc__)
  helpers.launch_datasets("scrnaseq/load_openproblems_op3", DATASETS, SHARED, "scrnaseq", args)
