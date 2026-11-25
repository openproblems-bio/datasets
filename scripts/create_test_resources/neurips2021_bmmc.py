# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "mudata==0.4.1", "pyyaml"]
# ///
"""Test resources: NeurIPS 2021 bone marrow.

Bone marrow mononuclear cells from healthy human donors, profiled with two multimodal
assays for the NeurIPS 2021 multimodal single-cell data integration challenge
(https://openproblems.bio/events/2021-09_neurips, Luecken et al., 2021):

- CITE-seq: gene expression (GEX) together with the abundance of cell surface proteins
  (ADT), measured with the 10x 3' kit and the BioLegend TotalSeq B Universal Human Panel.
- Multiome: gene expression together with chromatin accessibility (ATAC), measured with
  the 10x Multiome kit.

Samples were processed at four sites, and some donors were measured at several sites,
giving a nested batch structure.

Run from the root of the repository:

    uv run scripts/create_test_resources/neurips2021_bmmc.py
"""

import lamindb as ln
import mudata as md

import helpers

DATASETS = ["bmmc_cite", "bmmc_multiome"]

helpers.check_instance()
ln.track("aBQ5etnzMzCR", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/neurips2021_bmmc.py")

# load both datasets, subsampled to 600 cells, evenly across batches, and 1500 features of
# each modality (all of them for the surface proteins); this gives one raw MuData store per
# dataset, with modalities rna and prot or atac, and the cell metadata in the global obs
geo = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl"
common = {
  "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122",
  "dataset_reference": "luecken2021neurips",
}
raw_dir = helpers.run_workflow("multimodal/load_openproblems_neurips2021_bmmc", {
  "param_list": [
    {
      "id": "openproblems_neurips2021/bmmc_cite",
      "input": f"{geo}/GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz",
      "mod1": "GEX",
      "mod2": "ADT",
      "dataset_name": "OpenProblems NeurIPS2021 CITE-Seq",
      "dataset_summary": "Single-cell CITE-Seq (GEX+ADT) data collected from bone marrow mononuclear cells of 12 healthy human donors.",
      "dataset_description": (
        "Single-cell CITE-Seq data collected from bone marrow mononuclear cells of 12 healthy human donors "
        "using the 10X 3 prime Single-Cell Gene Expression kit with Feature Barcoding in combination with the "
        "BioLegend TotalSeq B Universal Human Panel v1.0. The dataset was generated to support Multimodal "
        "Single-Cell Data Integration Challenge at NeurIPS 2021. Samples were prepared using a standard protocol "
        "at four sites. The resulting data was then annotated to identify cell types and remove doublets. The "
        "dataset was designed with a nested batch layout such that some donor samples were measured at multiple "
        "sites with some donors measured at a single site."
      ),
      **common,
    },
    {
      "id": "openproblems_neurips2021/bmmc_multiome",
      "input": f"{geo}/GSE194122%5Fopenproblems%5Fneurips2021%5Fmultiome%5FBMMC%5Fprocessed%2Eh5ad%2Egz",
      "mod1": "GEX",
      "mod2": "ATAC",
      "dataset_name": "OpenProblems NeurIPS2021 Multiome",
      "dataset_summary": "Single-cell Multiome (GEX+ATAC) data collected from bone marrow mononuclear cells of 12 healthy human donors.",
      "dataset_description": (
        "Single-cell Multiome data collected from bone marrow mononuclear cells of 12 healthy human donors "
        "using the 10X Multiome Gene Expression and Chromatin Accessibility kit. The dataset was generated to "
        "support Multimodal Single-Cell Data Integration Challenge at NeurIPS 2021. Samples were prepared using "
        "a standard protocol at four sites. The resulting data was then annotated to identify cell types and "
        "remove doublets. The dataset was designed with a nested batch layout such that some donor samples were "
        "measured at multiple sites with some donors measured at a single site."
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
}, highmem="64.GB")

# normalise each modality in the way that suits it: gene expression to log CP10k, protein
# abundance with a centred log-ratio transformation (CLR), chromatin accessibility with
# TF-IDF; then compute a truncated SVD embedding and select highly variable features per
# modality, and combine the modalities into one MuData store again
output_dir = helpers.run_workflow("multimodal/process_multimodal", {
  "param_list": [
    {"id": f"openproblems_neurips2021/{dataset}", "input": str(raw_dir / "openproblems_neurips2021" / dataset / "raw.mudata.zarr")}
    for dataset in DATASETS
  ],
  "output_dataset": "$id/dataset.mudata.zarr",
  "output_state": "$id/state_process.yaml",
}, output_dir=raw_dir)

for dataset in DATASETS:
  for file in ["raw", "dataset"]:
    print(dataset, file, md.read_zarr(output_dir / "openproblems_neurips2021" / dataset / f"{file}.mudata.zarr"), sep="\n")

# store the raw and processed dataset as artifacts, in one collection
helpers.save_collection(
  output_dir / "openproblems_neurips2021",
  [f"{dataset}/{file}.mudata.zarr" for dataset in DATASETS for file in ["raw", "dataset"]],
  key="resources_test/datasets/openproblems_neurips2021",
  description="NeurIPS 2021 bone marrow CITE-seq and Multiome, subsampled to 600 cells and 1500 features per modality",
  modality="multimodal",
)

ln.finish()
