import anndata as ad
import mudata as md

## VIASH START
par = {
  "input_raw": "resources_test/common/openproblems_neurips2021/bmmc_cite/raw.mudata.zarr",
  "input_rna": "rna.anndata.zarr",
  "input_prot": "prot.anndata.zarr",
  "input_atac": None,
  "output": "dataset.mudata.zarr",
}
## VIASH END

print(">> Load data", flush=True)
raw = md.read_zarr(par["input_raw"])

mods = {}
for name in ["rna", "prot", "atac"]:
  path = par[f"input_{name}"]
  if name not in raw.mod:
    assert path is None, f"The raw dataset has no '{name}' modality, but --input_{name} was given"
    continue
  assert path is not None, f"The raw dataset has a '{name}' modality, but --input_{name} is missing"

  adata = ad.read_zarr(path)
  assert adata.obs_names.equals(raw.obs_names), f"Modality '{name}' has different cells than the raw dataset"

  # drop what split_mudata copied in from the global obs and uns, keep what processing added
  adata.obs = adata.obs.drop(columns=[col for col in raw.obs.columns if col in adata.obs])
  adata.uns = {
    **raw[name].uns,
    **{key: value for key, value in adata.uns.items() if key not in raw.uns},
  }
  mods[name] = adata

print(">> Combine the modalities", flush=True)
mdata = md.MuData(mods)
mdata.obs = raw.obs.copy()
mdata.uns = dict(raw.uns)
print(mdata, flush=True)

print(">> Writing data", flush=True)
mdata.write_zarr(par["output"])
