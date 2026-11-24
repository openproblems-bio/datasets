import subprocess
import sys

import anndata as ad
import mudata as md
import numpy as np

## VIASH START
meta = {
  "executable": "target/executable/utils/combine_mudata/combine_mudata",
  "resources_dir": "src/utils/split_mudata",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from test_data import make_mudata


def processed(mdata, name, extra_obs=True):
  """What a modality looks like after split_mudata, normalization, svd and hvg."""
  adata = mdata[name].copy()
  if extra_obs:
    adata.obs = adata.obs.join(mdata.obs)
  adata.uns = {**mdata.uns, **adata.uns, "normalization_id": f"norm_{name}"}
  adata.layers["normalized"] = adata.layers["counts"].astype(np.float64)
  adata.obs["size_factors"] = 1.0
  adata.obsm["X_svd"] = np.ones((adata.n_obs, 2))
  adata.var["hvg"] = True
  adata.var["hvg_score"] = 1.0
  path = f"{name}.zarr"
  adata.write_zarr(path)
  return path


print(">> Combine a CITE-seq-like store", flush=True)
raw = make_mudata("raw.mudata.zarr", mods=("rna", "prot"))
subprocess.run(
  [meta["executable"], "--input_raw", "raw.mudata.zarr", "--input_rna", processed(raw, "rna"), "--input_prot", processed(raw, "prot"), "--output", "dataset.mudata.zarr"],
  check=True,
)
mdata = md.read_zarr("dataset.mudata.zarr")

assert set(mdata.mod) == {"rna", "prot"}
assert list(mdata.obs.columns) == list(raw.obs.columns), f"global obs comes from the raw store, got {list(mdata.obs.columns)}"
assert all(":" not in col for col in mdata.obs.columns), "no prefixed columns in the global obs"
for name in ["rna", "prot"]:
  mod = mdata[name]
  assert list(mod.obs.columns) == [f"{name}_qc", "size_factors"], f"{name} obs: {list(mod.obs.columns)}"
  assert list(mod.var_names) == list(raw[name].var_names), f"{name} features are not renamed"
  assert mod.uns["normalization_id"] == f"norm_{name}"
  assert mod.uns["source"] == name, "the modality's own uns is kept"
  assert "dataset_id" not in mod.uns, "the copied global uns is dropped"
  for slot, key in [("layers", "normalized"), ("obsm", "X_svd"), ("var", "hvg"), ("var", "hvg_score")]:
    assert key in getattr(mod, slot), f"{name} misses {slot}['{key}']"
assert mdata.uns["dataset_id"] == "test"

print(">> Combine a store with only rna", flush=True)
raw = make_mudata("raw_rna.mudata.zarr", mods=("rna",))
subprocess.run(
  [meta["executable"], "--input_raw", "raw_rna.mudata.zarr", "--input_rna", processed(raw, "rna"), "--output", "dataset_rna.mudata.zarr"],
  check=True,
)
assert set(md.read_zarr("dataset_rna.mudata.zarr").mod) == {"rna"}

print(">> All tests passed successfully", flush=True)
