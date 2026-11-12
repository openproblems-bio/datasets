# Test resources

Each notebook in this directory creates one set of test resources: it runs a workflow of this repository on a small subset of a dataset and stores the outputs in the public LaminDB instance [`openproblems/data-dev`](https://lamin.ai/openproblems/data-dev), as one artifact per file grouped in a collection.

| Notebook | Collection |
|---|---|
| `pancreas.ipynb` | `resources_test/datasets/pancreas` |
| `cxg_immune_cell_atlas.ipynb` | `resources_test/datasets/cxg_immune_cell_atlas` |
| `cxg_mouse_pancreas_atlas.ipynb` | `resources_test/datasets/cxg_mouse_pancreas_atlas` |
| `mouse_brain_coronal.ipynb` | `resources_test/datasets/mouse_brain_coronal` |
| `neurips2021_bmmc.ipynb` | `resources_test/datasets/openproblems_neurips2021` |
| `scicar_cell_lines.ipynb` | `resources_test/datasets/scicar_cell_lines` |

`neurips2022_pbmc.sh` has not been converted yet, as its input data is not public.

## Running a notebook

You need Docker, Nextflow, the components built with `viash ns build --setup cachedbuild`, and write access to `openproblems/data-dev`. The `.lamin` directory in the root of this repository connects it to that instance. Run a notebook with:

```bash
uv run --with 'lamindb[aws]' --with nbconvert --with ipykernel \
  scripts/create_test_resources/run_notebook.sh scripts/create_test_resources/pancreas.ipynb
```

This tracks the run in LaminDB with the project `openproblems-bio/datasets`, links it to the notebook on GitHub if the notebook is committed and pushed, and saves the executed notebook as the report of the run.

The workflows are deterministic, so rerunning a notebook without changes returns the existing artifacts and collection. If the outputs change, the collection gets a new version; update its uid in `_viash.yaml` if this repository's own tests use it.

`helpers.py` holds the steps the notebooks share.
