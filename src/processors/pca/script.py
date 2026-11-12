import anndata as ad
import scanpy as sc

### VIASH START
par = {
  'input': 'resources_test/common/pancreas/dataset.anndata.zarr',
  'input_layer': 'log_cp10k',
  'output': 'dataset.zarr',
  'obsm_embedding': 'X_pca',
  'varm_loadings': 'pca_loadings',
  'uns_variance': 'pca_variance',
  'num_components': 25
}
### VIASH END

print(">> Load data", flush=True)
adata = ad.read_zarr(par['input'])

print(">> Look for layer", flush=True)
layer = adata.X if not par['input_layer'] else adata.layers[par['input_layer']]

print(">> Run PCA", flush=True)
X_pca, loadings, variance, variance_ratio = sc.tl.pca(
    layer, 
    n_comps=par["num_components"], 
    return_info=True
)

print(">> Storing output", flush=True)
adata.obsm[par["obsm_embedding"]] = X_pca
adata.varm[par["varm_loadings"]] = loadings.T
adata.uns[par["uns_variance"]] = {
    "variance": variance, 
    "variance_ratio": variance_ratio
}

print(">> Writing data", flush=True)
adata.write_zarr(par['output'])
