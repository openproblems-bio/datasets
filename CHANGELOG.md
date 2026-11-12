# datasets v0.2.0 (unreleased)

* MAJOR CHANGE: Bump Viash to 0.10.0-RC1 and rename `repository` to `package` in dependencies.

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
