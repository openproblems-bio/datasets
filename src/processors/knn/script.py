import anndata as ad
import scanpy as sc

### VIASH START
par = {
  'input': 'work/ca/0751ff85df6f9478cb7bda5a705cad/zebrafish.sqrt_cpm.pca.output.anndata.zarr',
  'input_layer': 'normalized',
  'output': 'dataset.zarr',
  'key_added': 'knn',
  'num_neighbors': 15,
  'seed': 0
}
### VIASH END

print(">> Load data", flush=True)
adata = ad.read_zarr(par['input'])

print(">> Run kNN", flush=True)
sc.pp.neighbors(
    adata,
    use_rep='X_pca',
    key_added=par['key_added'],
    n_neighbors=par['num_neighbors'],
    random_state=par['seed']
)

print(">> Writing data", flush=True)
adata.write_zarr(par['output'])
