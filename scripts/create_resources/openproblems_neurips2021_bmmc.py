# /// script
# requires-python = ">=3.12"
# dependencies = ["pyyaml"]
# ///
"""Load the NeurIPS 2021 multimodal datasets on de.NBI and publish them to LaminDB.

Run from the root of the repository:

    uv run scripts/create_resources/openproblems_neurips2021_bmmc.py [--ids ...] [--dry-run]
"""

import helpers

GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl"

DATASETS = [
  {
    "id": "openproblems_neurips2021/bmmc_cite",
    "input": f"{GEO}/GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz",
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
  },
  {
    "id": "openproblems_neurips2021/bmmc_multiome",
    "input": f"{GEO}/GSE194122%5Fopenproblems%5Fneurips2021%5Fmultiome%5FBMMC%5Fprocessed%2Eh5ad%2Egz",
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
  },
]

SHARED = {
  "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122",
  "dataset_reference": "luecken2021neurips",
  "output_raw": "$id/raw.mudata.zarr",
}

if __name__ == "__main__":
  args = helpers.parse_args(__doc__)
  helpers.launch_datasets("multimodal/load_openproblems_neurips2021_bmmc", DATASETS, SHARED, "multimodal", args)
