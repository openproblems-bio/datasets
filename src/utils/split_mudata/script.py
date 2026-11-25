import mudata as md

## VIASH START
par = {
  "input": "resources_test/common/openproblems_neurips2021/bmmc_cite/raw.mudata.zarr",
  "output_rna": "rna.anndata.zarr",
  "output_prot": "prot.anndata.zarr",
  "output_atac": "atac.anndata.zarr",
}
## VIASH END

print(">> Load data", flush=True)
mdata = md.read_zarr(par["input"])
print(mdata, flush=True)

unknown = set(mdata.mod) - {"rna", "prot", "atac"}
assert not unknown, f"Unknown modalities: {sorted(unknown)}"
assert "rna" in mdata.mod, "The dataset has no 'rna' modality"

for name, adata in mdata.mod.items():
  output = par[f"output_{name}"]
  if output is None:
    print(f">> Skipping modality '{name}', no output given", flush=True)
    continue

  # components for single-modality data read the cell metadata and dataset_* fields from the AnnData
  print(f">> Writing modality '{name}'", flush=True)
  adata = adata.copy()
  adata.obs = adata.obs.join(mdata.obs)
  adata.uns = {**adata.uns, **mdata.uns}
  adata.write_zarr(output)
