# Test resources

Each notebook or script in this directory creates one set of test resources: it runs a workflow of this repository on a small subset of a dataset and stores the outputs in the public LaminDB instance [`openproblems/data-dev`](https://lamin.ai/openproblems/data-dev), as one artifact per file grouped in a collection.

| Notebook or script | Collection |
|---|---|
| `pancreas.ipynb` | `resources_test/datasets/pancreas` |
| `cxg_immune_cell_atlas.ipynb` | `resources_test/datasets/cxg_immune_cell_atlas` |
| `cxg_mouse_pancreas_atlas.ipynb` | `resources_test/datasets/cxg_mouse_pancreas_atlas` |
| `mouse_brain_coronal.ipynb` | `resources_test/datasets/mouse_brain_coronal` |
| `neurips2021_bmmc.py` | `resources_test/datasets/openproblems_neurips2021` |
| `neurips2022_pbmc.py` | `resources_test/datasets/openproblems_neurips2022` |
| `scicar_cell_lines.ipynb` | `resources_test/datasets/scicar_cell_lines` |

The source files of `neurips2022_pbmc.py` are not public. They are stored in the private instance [`openproblems/private`](https://lamin.ai/openproblems/private), and `helpers.fetch_private()` records them in the parameters of the run rather than as its inputs, so that their records are not copied into the public instance.

`pancreas.ipynb` and `scicar_cell_lines.ipynb` use the deprecated openproblems v1 datasets. Their collections stay available, but new tests should use the CELLxGENE test resources instead.

## Running a script or notebook

You need Docker, Nextflow, the components built with `viash ns build --setup cachedbuild`, and write access to `openproblems/data-dev`. The `.lamin` directory in the root of this repository connects it to that instance, so run them from there. A script lists its dependencies, so `uv` installs them:

```bash
uv run scripts/create_test_resources/neurips2021_bmmc.py
```

Run a notebook with:

```bash
uv run --with 'lamindb[aws,bionty]==2.12.2' --with 'anndata==0.13.2' --with 'mudata==0.4.1' --with nbconvert --with ipykernel \
  scripts/create_test_resources/run_notebook.sh scripts/create_test_resources/pancreas.ipynb
```

This tracks the run in LaminDB with the project `openproblems-bio/datasets` and links it to the script or notebook on GitHub if it is committed and pushed. For a notebook, the executed notebook is saved as the report of the run.

The workflows are deterministic, so rerunning a script or notebook without changes returns the existing artifacts and collection. If the outputs change, the collection gets a new version; update its uid in `_viash.yaml` if this repository's own tests use it.

Each dataset artifact is also annotated with what the OpenProblems website shows about a dataset, with [`scripts/lamin/annotate_artifacts.py`](../lamin/annotate_artifacts.py).

`helpers.py` holds the steps they share.
