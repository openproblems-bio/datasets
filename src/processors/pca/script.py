import anndata as ad
import numpy as np
import scanpy as sc

### VIASH START
par = {
  'input': 'resources_test/common/cxg_mouse_pancreas_atlas/dataset.anndata.zarr',
  'input_layer': 'log_cp10k',
  'output': 'dataset.zarr',
  'obsm_embedding': 'X_pca',
  'varm_loadings': 'pca_loadings',
  'uns_variance': 'pca_variance',
  'num_components': 25,
  'input_var_features': 'hvg',
  'seed': 0
}
### VIASH END

print(">> Load data", flush=True)
adata = ad.read_zarr(par['input'])

print(">> Look for layer", flush=True)
layer = adata.X if not par['input_layer'] else adata.layers[par['input_layer']]

print(">> Select features", flush=True)
if par["input_var_features"]:
    feature_mask = adata.var[par["input_var_features"]].to_numpy(dtype=bool)
else:
    feature_mask = np.ones(adata.n_vars, dtype=bool)

print(">> Run PCA", flush=True)
X_pca, loadings, variance, variance_ratio = sc.tl.pca(
    layer[:, feature_mask],
    n_comps=par["num_components"],
    random_state=par["seed"],
    return_info=True
)

print(">> Storing output", flush=True)
adata.obsm[par["obsm_embedding"]] = X_pca
# features that were not selected get zero loadings
full_loadings = np.zeros((adata.n_vars, loadings.shape[0]), dtype=loadings.dtype)
full_loadings[feature_mask] = loadings.T
adata.varm[par["varm_loadings"]] = full_loadings
adata.uns[par["uns_variance"]] = {
    "variance": variance, 
    "variance_ratio": variance_ratio
}

print(">> Writing data", flush=True)
adata.write_zarr(par['output'])
