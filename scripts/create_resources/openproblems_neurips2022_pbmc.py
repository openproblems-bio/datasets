# /// script
# requires-python = ">=3.12"
# dependencies = ["pyyaml"]
# ///
"""Load the NeurIPS 2022 multimodal datasets on de.NBI and publish them to LaminDB.

The competition data, merged per modality, is not public: it is read from the private
instance openproblems/private, by artifact uid. The Multiome ATAC file is 16 GB, so that
dataset gets more memory and time.

Run from the root of the repository:

    uv run scripts/create_resources/openproblems_neurips2022_pbmc.py [--ids ...] [--dry-run]
"""

import helpers

PRIVATE = "lamin://openproblems/private/artifact"

DATASETS = [
  {
    "id": "openproblems_neurips2022/pbmc_cite",
    "input_mod1": f"{PRIVATE}/0ReyMrCccBLuMgqT0000",  # neurips2022/cite_rna_merged.h5ad
    "input_mod2": f"{PRIVATE}/Wg9Ggwp6rBSnxVSK0000",  # neurips2022/cite_prot_merged.h5ad
    "mod1": "GEX",
    "mod2": "ADT",
    "dataset_name": "OpenProblems NeurIPS2022 CITE-Seq",
    "dataset_summary": "Single-cell CITE-Seq (GEX+ADT) data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells from 4 healthy human donors, over 4 days.",
    "dataset_description": (
      "Single-cell CITE-Seq data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor "
      "cells (HSPCs) from 4 healthy human donors, measured on days 2, 3, 4 and 7 of a 10-day time course. "
      "Gene expression was measured together with the abundance of 140 cell surface proteins. The dataset "
      "was generated for the Multimodal Single-Cell Integration competition at NeurIPS 2022, in which the "
      "abundance of the proteins had to be predicted from gene expression on a later day. The cells were "
      "annotated with seven progenitor cell types."
    ),
  },
  {
    "id": "openproblems_neurips2022/pbmc_multiome",
    "input_mod1": f"{PRIVATE}/iYN6PKlHwv6hiM1B0000",  # neurips2022/multiome_rna_merged.h5ad
    "input_mod2": f"{PRIVATE}/GZUenxpFeaPXjwjX0000",  # neurips2022/multiome_atac_merged.h5ad
    "mod1": "GEX",
    "mod2": "ATAC",
    "dataset_name": "OpenProblems NeurIPS2022 Multiome",
    "dataset_summary": "Single-cell Multiome (GEX+ATAC) data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells from 4 healthy human donors, over 5 days.",
    "dataset_description": (
      "Single-cell Multiome data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor "
      "cells (HSPCs) from 4 healthy human donors, measured on days 2, 3, 4, 7 and 10 of a 10-day time "
      "course with the 10X Multiome Gene Expression and Chromatin Accessibility kit. The dataset was "
      "generated for the Multimodal Single-Cell Integration competition at NeurIPS 2022, in which gene "
      "expression had to be predicted from chromatin accessibility on a later day. The cells were annotated "
      "with seven progenitor cell types; the cell types of the cells in the private test set were not "
      "released."
    ),
    "large": True,
  },
]

SHARED = {
  "dataset_url": "https://www.kaggle.com/competitions/open-problems-multimodal/data",
  "dataset_reference": "lance2024predicting",
  "output_raw": "$id/raw.mudata.zarr",
  "output_state": "$id/load_state.yaml",
}

if __name__ == "__main__":
  args = helpers.parse_args(__doc__)
  helpers.launch_datasets("multimodal/load_openproblems_neurips2022_pbmc", DATASETS, SHARED, "multimodal", args)
