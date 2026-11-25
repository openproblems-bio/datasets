# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "mudata==0.4.1", "pyyaml"]
# ///
"""Test resources: NeurIPS 2022 CD34+ HSPCs.

Mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells (HSPCs) from four
healthy human donors, profiled over a 10-day time course for the NeurIPS 2022 multimodal
single-cell integration competition on Kaggle
(https://www.kaggle.com/competitions/open-problems-multimodal):

- CITE-seq: gene expression (GEX) together with the abundance of 140 cell surface proteins
  (ADT), on days 2, 3, 4 and 7.
- Multiome: gene expression together with chromatin accessibility (ATAC), on days 2, 3, 4,
  7 and 10.

The cells are annotated with seven progenitor cell types. The cell types of the Multiome
cells in the private test set of the competition were not released, and are labelled
`hidden`.

The source files are stored in the private LaminDB instance `openproblems/private`, so you
need access to it to run this script. The test resources themselves are public.

Run from the root of the repository:

    uv run scripts/create_test_resources/neurips2022_pbmc.py
"""

import lamindb as ln
import mudata as md

import helpers

DATASETS = ["pbmc_cite", "pbmc_multiome"]

helpers.check_instance()
ln.track("fTZVtuOVzHtH", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/neurips2022_pbmc.py")

# the training and test data of the competition, merged per modality
inputs = {
  file: helpers.fetch_private(f"neurips2022/{file}.h5ad")
  for file in ["cite_rna_merged", "cite_prot_merged", "multiome_rna_merged", "multiome_atac_merged"]
}

# load both datasets, subsampled to 600 cells, evenly across batches, and 1500 features of
# each modality (all of them for the surface proteins); this gives one raw MuData store per
# dataset, with modalities rna and prot or atac, and the cell metadata in the global obs
common = {
  "dataset_url": "https://www.kaggle.com/competitions/open-problems-multimodal/data",
  "dataset_reference": "lance2024predicting",
}
raw_dir = helpers.run_workflow("multimodal/load_openproblems_neurips2022_pbmc", {
  "param_list": [
    {
      "id": "openproblems_neurips2022/pbmc_cite",
      "input_mod1": str(inputs["cite_rna_merged"]),
      "input_mod2": str(inputs["cite_prot_merged"]),
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
      **common,
    },
    {
      "id": "openproblems_neurips2022/pbmc_multiome",
      "input_mod1": str(inputs["multiome_rna_merged"]),
      "input_mod2": str(inputs["multiome_atac_merged"]),
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
      **common,
    },
  ],
  "do_subsample": True,
  "even": True,
  "n_obs": 600,
  "n_vars": 1500,
  "seed": 123,
  "output_raw": "$id/raw.mudata.zarr",
  "output_state": "$id/state.yaml",
}, highmem="100.GB")

# normalise each modality in the way that suits it: gene expression to log CP10k, protein
# abundance with a centred log-ratio transformation (CLR), chromatin accessibility with
# TF-IDF; then compute a truncated SVD embedding and select highly variable features per
# modality, and combine the modalities into one MuData store again
output_dir = helpers.run_workflow("multimodal/process_multimodal", {
  "param_list": [
    {"id": f"openproblems_neurips2022/{dataset}", "input": str(raw_dir / "openproblems_neurips2022" / dataset / "raw.mudata.zarr")}
    for dataset in DATASETS
  ],
  "output_dataset": "$id/dataset.mudata.zarr",
  "output_state": "$id/state_process.yaml",
}, output_dir=raw_dir)

for dataset in DATASETS:
  for file in ["raw", "dataset"]:
    print(dataset, file, md.read_zarr(output_dir / "openproblems_neurips2022" / dataset / f"{file}.mudata.zarr"), sep="\n")

# store the raw and processed dataset as artifacts, in one collection
helpers.save_collection(
  output_dir / "openproblems_neurips2022",
  [f"{dataset}/{file}.mudata.zarr" for dataset in DATASETS for file in ["raw", "dataset"]],
  key="resources_test/datasets/openproblems_neurips2022",
  description="NeurIPS 2022 CD34+ HSPC CITE-seq and Multiome, subsampled to 600 cells and 1500 features per modality",
  modality="multimodal",
)

ln.finish()
