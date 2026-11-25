# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "pyyaml"]
# ///
"""Test resources: Immune Cell Atlas.

Immune cells from 16 tissues of 12 adult donors, profiled with scRNA-seq (Domínguez Conde
et al., 2022, https://doi.org/10.1126/science.abl5197). The atlas was used to train
CellTypist and covers about 360,000 cells with finely annotated immune cell types.

Run from the root of the repository:

    uv run scripts/create_test_resources/cxg_immune_cell_atlas.py
"""

import anndata as ad
import lamindb as ln

import helpers

helpers.check_instance()
ln.track("RXx1ZGPcLsGt", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/cxg_immune_cell_atlas.py")

# take the cells of three donors (A29, A31 and A35) from the CELLxGENE Census release of
# 2023-07-25, use the donor as batch, and subsample to 600 cells and 1500 genes
output_dir = helpers.run_workflow("scrnaseq/load_cellxgene_census", {
  "id": "cxg_immune_cell_atlas",
  "species": "homo_sapiens",
  "census_version": "2023-07-25",
  "obs_value_filter": "dataset_id == '1b9d8702-5af8-4142-85ed-020eb06ec4f6' and donor_id in ['A29', 'A31', 'A35']",
  "obs_batch": ["donor_id"],
  "dataset_name": "Immune Cell Atlas",
  "dataset_summary": "Cross-tissue immune cell analysis reveals tissue-specific features in humans",
  "dataset_description": (
    "Despite their crucial role in health and disease, our knowledge of immune cells within human tissues "
    "remains limited. We surveyed the immune compartment of 16 tissues from 12 adult donors by single-cell "
    "RNA sequencing and VDJ sequencing generating a dataset of ~360,000 cells. To systematically resolve "
    "immune cell heterogeneity across tissues, we developed CellTypist, a machine learning tool for rapid "
    "and precise cell type annotation. Using this approach, combined with detailed curation, we determined "
    "the tissue distribution of finely phenotyped immune cell types, revealing hitherto unappreciated "
    "tissue-specific features and clonal architecture of T and B cells. Our multitissue approach lays the "
    "foundation for identifying highly resolved immune cell types by leveraging a common reference dataset, "
    "tissue-integrated expression analysis, and antigen receptor sequencing."
  ),
  "dataset_url": "https://cellxgene.cziscience.com/collections/62ef75e4-cbea-454e-a0ce-998ec40223d3",
  "dataset_reference": "dominguez2022crosstissue",
  "do_subsample": True,
  "n_obs": 600,
  "n_vars": 1500,
  "seed": 123,
  "output_raw": "$id/raw.anndata.zarr",
  "output_state": "$id/load_state.yaml",
})

# normalise the counts to log CP10k, select the 1000 most highly variable genes, compute a
# 50-dimensional PCA on those genes and build a 15-nearest-neighbour graph on the PCA
output_dir = helpers.run_workflow("scrnaseq/process_scrnaseq", {
  "id": "cxg_immune_cell_atlas",
  "input": str(output_dir / "cxg_immune_cell_atlas" / "raw.anndata.zarr"),
  "output_dataset": "$id/dataset.anndata.zarr",
  "output_state": "$id/process_state.yaml",
}, output_dir=output_dir)

print(ad.read_zarr(output_dir / "cxg_immune_cell_atlas" / "dataset.anndata.zarr"))

# store the processed dataset as an artifact, in a collection
helpers.save_collection(
  output_dir / "cxg_immune_cell_atlas",
  ["dataset.anndata.zarr"],
  key="resources_test/datasets/cxg_immune_cell_atlas",
  modality="scrnaseq",
  description="Immune Cell Atlas (CELLxGENE), three donors, subsampled to 600 cells and 1500 genes",
)

ln.finish()
