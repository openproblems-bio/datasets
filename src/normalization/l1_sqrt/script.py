import anndata as ad
import numpy as np
import scipy.sparse as sp

## VIASH START
par = {
    'input': "output_train.anndata.zarr",
    'output': "output_norm.zarr"
}
meta = {
    'name': "l1_sqrt"
}
## VIASH END

print("Load data", flush=True)
adata = ad.read_zarr(par['input'])

print("Normalize data", flush=True)
# sqrt of the counts
sqrt_data = sp.csr_matrix(adata.layers['counts'], copy=True)
sqrt_data.data = np.sqrt(sqrt_data.data)

# l1 normalize each cell, leaving empty cells at zero
libsize = np.asarray(sqrt_data.sum(axis=1)).ravel()
l1_sqrt = sp.diags(1 / np.where(libsize == 0, 1, libsize)) @ sqrt_data
l1_sqrt = l1_sqrt.tocsr()

print("Store output in adata", flush=True)
adata.layers[par["layer_output"]] = l1_sqrt
adata.uns["normalization_id"] = par["normalization_id"] or meta['name']

print("Write data", flush=True)
adata.write_zarr(par['output'])
