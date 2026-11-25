# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "pyyaml"]
# ///
"""Test resources: mouse pancreatic islet atlas.

A mouse pancreatic islet atlas of over 300,000 cells from nine datasets, spanning sexes,
ages from embryonic to aged, chemical stress and models of type 1 and type 2 diabetes
(Hrovatin et al., 2023, https://doi.org/10.1038/s42255-023-00876-x).

The output of every processing step is kept, so that each component can be tested on the
input it gets in the workflow. The component tests of this repository use these.

Run from the root of the repository:

    uv run scripts/create_test_resources/cxg_mouse_pancreas_atlas.py
"""

import anndata as ad
import lamindb as ln

import helpers

STEPS = ["raw", "normalized", "hvg", "pca", "dataset"]

helpers.check_instance()
ln.track("xLOCoaLOtD9h", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/cxg_mouse_pancreas_atlas.py")

# take the cells of three two-year-old mice from the CELLxGENE Census release of
# 2025-11-08, use the donor as batch, and subsample to 600 cells and 1500 genes
donors = [
  "mouse_pancreatic_islet_atlas_Hrovatin__Fltp_2y__MUC13974",
  "mouse_pancreatic_islet_atlas_Hrovatin__Fltp_2y__MUC13975",
  "mouse_pancreatic_islet_atlas_Hrovatin__Fltp_2y__MUC13976",
]
output_dir = helpers.run_workflow("scrnaseq/load_cellxgene_census", {
  "id": "cxg_mouse_pancreas_atlas",
  "species": "mus_musculus",
  "census_version": "2025-11-08",
  "obs_value_filter": f"dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8' and donor_id in {donors}",
  "obs_batch": ["donor_id"],
  "dataset_name": "Mouse Pancreatic Islet Atlas",
  "dataset_summary": "Mouse pancreatic islet scRNA-seq atlas across sexes, ages, and stress conditions including diabetes",
  "dataset_description": (
    "To better understand pancreatic β-cell heterogeneity we generated a mouse pancreatic islet atlas "
    "capturing a wide range of biological conditions. The atlas contains scRNA-seq datasets of over 300,000 "
    "mouse pancreatic islet cells, of which more than 100,000 are β-cells, from nine datasets with 56 "
    "samples, including two previously unpublished datasets. The samples vary in sex, age (ranging from "
    "embryonic to aged), chemical stress, and disease status (including T1D NOD model development and two "
    "T2D models, mSTZ and db/db) together with different diabetes treatments. Additional information about "
    "data fields is available in anndata uns field 'field_descriptions' and on "
    "https://github.com/theislab/mm_pancreas_atlas_rep/blob/main/resources/cellxgene.md."
  ),
  "dataset_url": "https://cellxgene.cziscience.com/collections/296237e2-393d-4e31-b590-b03f74ac5070",
  "dataset_reference": "hrovatin2023delineating",
  "do_subsample": True,
  "n_obs": 600,
  "n_vars": 1500,
  "seed": 123,
  "output_raw": "$id/raw.anndata.zarr",
  "output_state": "$id/load_state.yaml",
})

# normalise the counts to log CP10k, select the 1000 most highly variable genes, compute a
# 50-dimensional PCA on those genes and build a 15-nearest-neighbour graph on the PCA,
# keeping the output of every step
output_dir = helpers.run_workflow("scrnaseq/process_scrnaseq", {
  "id": "cxg_mouse_pancreas_atlas",
  "input": str(output_dir / "cxg_mouse_pancreas_atlas" / "raw.anndata.zarr"),
  **{f"output_{step}": f"$id/{step}.anndata.zarr" for step in STEPS[1:]},
  "output_state": "$id/process_state.yaml",
}, output_dir=output_dir)

print(ad.read_zarr(output_dir / "cxg_mouse_pancreas_atlas" / "dataset.anndata.zarr"))

# store the output of every step as an artifact, in one collection
helpers.save_collection(
  output_dir / "cxg_mouse_pancreas_atlas",
  [f"{step}.anndata.zarr" for step in STEPS],
  key="resources_test/datasets/cxg_mouse_pancreas_atlas",
  modality="scrnaseq",
  description="Mouse pancreatic islet atlas (CELLxGENE), three donors, subsampled to 600 cells and 1500 genes",
)

ln.finish()
