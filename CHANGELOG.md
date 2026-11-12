# datasets v0.2.0 (unreleased)

* MINOR CHANGE: Bump Viash to 0.9.7.

* MINOR CHANGE: Failed Nextflow processes are no longer ignored, so a failing dataset fails the run.

* MAJOR CHANGE: All loaders, components and workflows write AnnData zarr stores instead of h5ad files, and the file formats in `src/api` are now `anndata_zarr`. Every Docker image sets `ANNDATA_ZARR_WRITE_FORMAT=3`. Multimodal datasets are still two AnnData stores, spatial datasets one.

* MAJOR CHANGE: Remove the `extract_uns_metadata` step and `--output_meta*` outputs from the workflows, as well as the `extract_dataset_info` and `extract_dataset_meta` workflows. Dataset metadata will be stored in LaminDB.

* MINOR CHANGE: Remove the h5ad-only `--output_compression` argument.

* MINOR CHANGE: The openproblems v1 loaders install openproblems v0.8.0 without its pinned 2022-era dependencies and run on `openproblems/base_python:1`, so they write the same zarr format as the other loaders.

* BUG FIX: `processors/subsample` draws with a seeded numpy generator; `--seed` had no effect before.

* BUG FIX: `processors/pca` runs on the features selected by `--input_var_features` (HVGs by default), which it previously ignored. `processors/pca` and `processors/knn` gain a `--seed` argument.

* MINOR CHANGE: Update the `common` submodule, which adds `labels_denbi.config`.

* MINOR CHANGE: Drop the custom `multiple_sep` from the CELLxGENE census loaders and workflow; pass multiple values with `;`.

* MINOR CHANGE: Remove the unused `scripts/run_benchmark` scripts.

* DOCUMENTATION: Cite the published OpenProblems paper.

* BUG FIX: `normalization/l1_sqrt` computes the normalization with scipy instead of the unmaintained `scprep`, which no longer installs.

* BUG FIX: `processors/svd` uses `anndata.read_h5ad()` instead of the removed `anndata.read()`.

* BUG FIX: The openproblems v1 loaders build on Python 3.10 and download without `scprep`'s spoofed user agent, which figshare and GEO refuse.

# datasets v0.1.0

* Port dataset loaders and processors from `openproblems-bio/openproblems`.

* Port OP3 data loader (PR #2).

* Apply fixes to scripts and components (PR #3).
