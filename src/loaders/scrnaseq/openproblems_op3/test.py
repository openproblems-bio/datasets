import sys

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

## VIASH START
meta = {
  "executable": "target/executable/loaders/scrnaseq/openproblems_op3/openproblems_op3",
}
## VIASH END

DONORS = ["Donor 1", "Donor 2", "Donor 3"]
CELL_TYPES = ["B cells", "T cells", "Myeloid cells", "NK cells"]
COMPOUNDS = ["Dimethyl Sulfoxide", "Dabrafenib"]


@pytest.fixture
def source_h5ad(tmp_path):
  """A small h5ad shaped like GSE279945_sc_counts_processed.h5ad, analysis results included."""
  rng = np.random.default_rng(1)
  n_obs, n_vars = 60, 20
  obs = pd.DataFrame(
    {
      "sm_name": rng.choice(COMPOUNDS, n_obs),
      "donor_id": rng.choice(DONORS, n_obs),
      "cell_type": rng.choice(CELL_TYPES, n_obs),
      "plate_name": "plate_1",
      "leiden_res1": rng.choice(["0", "1"], n_obs),
    },
    index=[f"cell_{i}" for i in range(n_obs)],
  )
  obs["cell_type_orig"] = obs["cell_type"]
  obs["split"] = np.where(obs["sm_name"] == "Dimethyl Sulfoxide", "control", "train")
  var = pd.DataFrame(
    {"highly_variable": rng.random(n_vars) > 0.5, "means": rng.random(n_vars), "n_cells": rng.integers(1, 60, n_vars)},
    index=[f"GENE{i}" for i in range(n_vars)],
  )
  counts = sp.csr_matrix(rng.poisson(1, (n_obs, n_vars)).astype(np.float32))
  adata = ad.AnnData(X=counts, obs=obs, var=var)
  adata.raw = adata
  adata.obsm["X_pca"] = rng.random((n_obs, 5))
  adata.obsm["HTO_clr"] = pd.DataFrame(rng.random((n_obs, 2)), index=obs.index, columns=["HTO1", "HTO2"])
  adata.obsp["connectivities"] = sp.random(n_obs, n_obs, density=0.1, format="csr")
  adata.uns["log1p"] = {"base": None}
  adata.uns["leiden_res1_colors"] = ["#000000", "#ffffff"]
  adata.uns["neighbors"] = {"params": {"n_neighbors": 15}}
  path = tmp_path / "source.h5ad"
  adata.write_h5ad(path)
  return path


def test_loader(run_component, source_h5ad, tmp_path):
  output = tmp_path / "output.zarr"
  run_component([
    "--input", str(source_h5ad),
    "--donor_id", "Donor 1",
    "--dataset_id", "test_op3",
    "--dataset_name", "OP3 test",
    "--dataset_summary", "Test summary",
    "--dataset_description", "Test description",
    "--dataset_url", "https://example.org",
    "--dataset_reference", "GSE279945",
    "--output", str(output),
  ])

  adata = ad.read_zarr(output)

  assert adata.n_obs > 0
  assert set(adata.obs["donor_id"]) == {"Donor 1"}
  for column in ["cell_type_ontology_term_id", "sex_ontology_term_id", "sm_name", "split", "control"]:
    assert column in adata.obs

  # the counts and the cell and gene metadata, but not the analysis results of the source
  assert adata.X is None
  assert list(adata.layers) == ["counts"]
  assert adata.layers["counts"].sum() > 0
  assert list(adata.var.columns) == ["feature_name"]
  assert list(adata.var["feature_name"]) == list(adata.var_names)
  assert "leiden_res1" not in adata.obs
  assert adata.raw is None
  assert list(adata.obsm) == []
  assert list(adata.obsp) == []
  assert sorted(adata.uns) == sorted([
    "dataset_id", "dataset_name", "dataset_summary", "dataset_description",
    "dataset_url", "dataset_reference", "dataset_version",
  ])
  assert adata.uns["dataset_id"] == "test_op3"


if __name__ == "__main__":
  sys.exit(pytest.main([__file__]))
