# Test resources

Each script in this directory creates one set of test resources: it runs a workflow of this repository on a small subset of a dataset and stores the outputs in the public LaminDB instance [`openproblems/data-dev`](https://lamin.ai/openproblems/data-dev), as one artifact per file grouped in a collection.

| Script | Collection |
|---|---|
| `pancreas.py` | `resources_test/datasets/pancreas` |
| `cxg_immune_cell_atlas.py` | `resources_test/datasets/cxg_immune_cell_atlas` |
| `cxg_mouse_pancreas_atlas.py` | `resources_test/datasets/cxg_mouse_pancreas_atlas` |
| `mouse_brain_coronal.py` | `resources_test/datasets/mouse_brain_coronal` |
| `neurips2021_bmmc.py` | `resources_test/datasets/openproblems_neurips2021` |
| `neurips2022_pbmc.py` | `resources_test/datasets/openproblems_neurips2022` |
| `scicar_cell_lines.py` | `resources_test/datasets/scicar_cell_lines` |

The source files of `neurips2022_pbmc.py` are not public. They are stored in the private instance [`openproblems/private`](https://lamin.ai/openproblems/private), and `helpers.fetch_private()` records them in the parameters of the run rather than as its inputs, so that their records are not copied into the public instance.

`pancreas.py` and `scicar_cell_lines.py` use the deprecated openproblems v1 datasets. Their collections stay available, but new tests should use the CELLxGENE and NeurIPS test resources instead.

## Running a script

You need Docker, Nextflow, the components built with `viash ns build --setup cachedbuild`, and write access to `openproblems/data-dev`. The `.lamin` directory in the root of this repository connects it to that instance, so run the scripts from there. Each script lists its dependencies, so `uv` installs them:

```bash
uv run scripts/create_test_resources/cxg_mouse_pancreas_atlas.py
```

This tracks the run in LaminDB with the project `openproblems-bio/datasets` and links it to the script on GitHub if the script is committed and pushed.

The workflows are deterministic, so rerunning a script without changes returns the existing artifacts and collection. If the outputs change, the collection gets a new version; update its uid in `_viash.yaml` if this repository's own tests use it.

Each dataset artifact is also annotated with what the OpenProblems website shows about a dataset, with [`scripts/lamin/annotate_artifacts.py`](../lamin/annotate_artifacts.py).

`helpers.py` holds the steps the scripts share.
