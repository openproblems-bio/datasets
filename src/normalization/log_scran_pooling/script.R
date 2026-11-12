cat(">> Loading dependencies\n")
library(anndata, warn.conflicts = FALSE)
requireNamespace("scran", quietly = TRUE)
requireNamespace("BiocParallel", quietly = TRUE)
library(Matrix, warn.conflicts = FALSE)

## VIASH START
par <- list(
  input = "resources_test/common/pancreas/raw.anndata.zarr",
  output = "output.scran.zarr",
  layer_output = "log_scran_pooling",
  obs_size_factors = "size_factors_log_scran_pooling"
)
## VIASH END

cat(">> Load data\n")
# anndata::read_zarr() calls itself in R anndata 0.8.0, read through python instead
adata <- reticulate::import("anndata")$read_zarr(par$input)
counts <- as(t(adata$layers[["counts"]]), "CsparseMatrix")

cat(">> Normalizing data\n")
size_factors <- scran::calculateSumFactors(
  counts,
  min.mean = 0.1,
  BPPARAM = BiocParallel::MulticoreParam()
)
lognorm <- log1p(sweep(adata$layers[["counts"]], 1, size_factors, "*"))

cat(">> Storing in anndata\n")
adata$obs[[par$obs_size_factors]] <- size_factors
adata$layers[[par$layer_output]] <- lognorm
norm_id <- par[["normalization_id"]]
if (is.null(norm_id)) {
  norm_id <- meta[["name"]]
}
adata$uns[["normalization_id"]] <- norm_id

cat(">> Writing to file\n")
zzz <- adata$write_zarr(par$output)
