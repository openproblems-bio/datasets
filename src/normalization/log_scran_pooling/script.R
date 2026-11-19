cat(">> Loading dependencies\n")
requireNamespace("reticulate", quietly = TRUE)
requireNamespace("scran", quietly = TRUE)
requireNamespace("BiocParallel", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/common/cxg_mouse_pancreas_atlas/raw.anndata.zarr",
  output = "output.scran.zarr",
  layer_output = "log_scran_pooling",
  obs_size_factors = "size_factors_log_scran_pooling"
)
## VIASH END

cat(">> Load data\n")
# read through python anndata, since anndataR can't read sharded zarr yet
ad <- reticulate::import("anndata", convert = FALSE)
adata <- ad$read_zarr(par$input)
counts <- as(t(reticulate::py_to_r(adata$layers[["counts"]])), "CsparseMatrix")

cat(">> Normalizing data\n")
size_factors <- scran::calculateSumFactors(
  counts,
  min.mean = 0.1,
  BPPARAM = BiocParallel::MulticoreParam()
)
lognorm <- log1p(Matrix::Diagonal(x = size_factors) %*% Matrix::t(counts))

cat(">> Storing in anndata\n")
norm_id <- par[["normalization_id"]]
if (is.null(norm_id)) {
  norm_id <- meta[["name"]]
}
reticulate::py_set_item(adata$obs, par$obs_size_factors, reticulate::r_to_py(size_factors))
reticulate::py_set_item(adata$layers, par$layer_output, reticulate::r_to_py(as(lognorm, "CsparseMatrix")))
reticulate::py_set_item(adata$uns, "normalization_id", norm_id)

cat(">> Writing to file\n")
adata$write_zarr(par$output)
