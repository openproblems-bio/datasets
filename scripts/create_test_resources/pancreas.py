# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "pyyaml"]
# ///
"""Test resources: human pancreas.

Human pancreatic islet cells from six scRNA-seq studies, each using a different technology
(CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1 and SMARTer-seq), as collected for the
scIB benchmark (Luecken et al., 2022, https://doi.org/10.1038/s41592-021-01336-8). The
technologies leave strong batch effects, which makes this a classic test case for data
integration.

This uses the deprecated openproblems v1 datasets. The collection stays available, but new
tests should use the CELLxGENE test resources instead.

Run from the root of the repository:

    uv run scripts/create_test_resources/pancreas.py
"""

import urllib.request

import anndata as ad
import lamindb as ln

import helpers

STEPS = ["raw", "normalized", "hvg", "pca", "dataset"]

helpers.check_instance()
ln.track("GKVzEhqMMjfX", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/pancreas.py")

# the S-phase and G2/M-phase marker genes of Tirosh et al. (2016), as used by scIB, are
# always kept, so that cell cycle effects can still be scored
scib_url = "https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources"
cell_cycle_genes = []
for name in ["s_genes_tirosh_hm.txt", "g2m_genes_tirosh_hm.txt"]:
  with urllib.request.urlopen(f"{scib_url}/{name}") as response:
    cell_cycle_genes += response.read().decode().split()

# keep the acinar and beta cells measured with three of the technologies (CEL-seq, inDrop
# and SMARTer-seq), subsample to 600 cells and 1500 genes, normalise to log CP10k, select
# the 1000 most highly variable genes, compute a 50-dimensional PCA on those genes and
# build a 15-nearest-neighbour graph on the PCA, keeping the output of every step
output_dir = helpers.run_workflow("scrnaseq/process_openproblems_v1", {
  "id": "pancreas",
  "input_id": "pancreas",
  "obs_cell_type": "celltype",
  "obs_batch": "tech",
  "var_feature_name": "index",
  "layer_counts": "counts",
  "dataset_name": "Human pancreas",
  "dataset_url": "https://theislab.github.io/scib-reproducibility/dataset_pancreas.html",
  "dataset_reference": "luecken2022benchmarking",
  "dataset_summary": "Human pancreas cells dataset from the scIB benchmarks",
  "dataset_description": "Human pancreatic islet scRNA-seq data from 6 datasets across technologies (CEL-seq, CEL-seq2, Smart-seq2, inDrop, Fluidigm C1, and SMARTER-seq).",
  "dataset_organism": "homo_sapiens",
  "keep_cell_type_categories": ["acinar", "beta"],
  "keep_batch_categories": ["celseq", "inDrop4", "smarter"],
  "keep_features": cell_cycle_genes,
  "do_subsample": True,
  "n_obs": 600,
  "n_vars": 1500,
  "seed": 123,
  "normalization_methods": ["log_cp10k"],
  **{f"output_{step}": f"$id/{step}.anndata.zarr" for step in STEPS},
  "output_state": "$id/state.yaml",
})

print(ad.read_zarr(output_dir / "pancreas" / "dataset.anndata.zarr"))

# store the output of every step as an artifact, in one collection
helpers.save_collection(
  output_dir / "pancreas",
  [f"{step}.anndata.zarr" for step in STEPS],
  key="resources_test/datasets/pancreas",
  description="Human pancreas (scIB), subsampled to 600 cells and 1500 genes",
)

ln.finish()
