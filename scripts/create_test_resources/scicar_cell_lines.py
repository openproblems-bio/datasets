# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "pyyaml"]
# ///
"""Test resources: sci-CAR cell lines.

sci-CAR jointly measures the transcriptome and chromatin accessibility of the same cells
using combinatorial indexing (Cao et al., 2018, https://doi.org/10.1126/science.aau0730).
This dataset contains 4,825 cells from three cell lines (human HEK293T and A549 cells, and
mouse NIH/3T3 cells), collected 0, 1 and 3 hours after dexamethasone treatment.

This uses the deprecated openproblems v1 datasets, and stores the two modalities as two
AnnData stores. The collection stays available, but new tests should use the NeurIPS
multimodal test resources instead.

Run from the root of the repository:

    uv run scripts/create_test_resources/scicar_cell_lines.py
"""

import anndata as ad
import lamindb as ln

import helpers

MODALITIES = ["mod1", "mod2"]

helpers.check_instance()
ln.track("Aq0TwdUdIrHt", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/scicar_cell_lines.py")

# subsample to 600 cells and 1500 features per modality, normalise the gene expression to
# log CP10k and the chromatin accessibility with TF-IDF, and compute a truncated SVD
# embedding and select highly variable features for each modality
output_dir = helpers.run_workflow("multimodal/process_openproblems_v1_multimodal", {
  "id": "scicar_cell_lines",
  "input_id": "scicar_cell_lines",
  "obs_tissue": "source",
  "layer_counts": "counts",
  "obs_cell_type": "cell_name",
  "var_feature_id": "index",
  "var_feature_name": "gene_short_name",
  "dataset_name": "sci-CAR cell lines",
  "dataset_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117089",
  "dataset_reference": "cao2018joint",
  "dataset_summary": "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the accessibility of cellular chromatin in the same cells",
  "dataset_description": (
    "sciCAR is a combinatorial indexing-based assay that jointly measures cellular transcriptomes and the "
    "accessibility of cellular chromatin in the same cells. Here, we use two sciCAR datasets that were "
    "obtained from the same study. The first dataset contains 4,825 cells from three cell lines (HEK293T "
    "cells, NIH/3T3 cells, and A549 cells) at multiple timepoints (0, 1 hour, 3 hours) after dexamethasone "
    "treatment. The second dataset contains 11,233 cells from wild-type adult mouse kidney."
  ),
  "dataset_organism": ["homo_sapiens", "mus_musculus"],
  "mod1": "GEX",
  "mod2": "ATAC",
  "do_subsample": True,
  "n_obs": 600,
  "n_vars": 1500,
  "seed": 123,
  "normalization_methods": ["log_cp10k"],
  "output_mod1": "$id/dataset_mod1.anndata.zarr",
  "output_mod2": "$id/dataset_mod2.anndata.zarr",
  "output_state": "$id/state.yaml",
})

for mod in MODALITIES:
  print(mod, ad.read_zarr(output_dir / "scicar_cell_lines" / f"dataset_{mod}.anndata.zarr"), sep="\n")

# store both modalities as artifacts, in one collection
helpers.save_collection(
  output_dir / "scicar_cell_lines",
  [f"dataset_{mod}.anndata.zarr" for mod in MODALITIES],
  key="resources_test/datasets/scicar_cell_lines",
  description="sci-CAR cell lines, subsampled to 600 cells and 1500 features per modality",
)

ln.finish()
