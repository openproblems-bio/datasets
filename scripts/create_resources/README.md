# Resources

The scripts in this directory build the full OpenProblems datasets on the de.NBI cluster and publish them to LaminDB with [nf-lamin](https://github.com/laminlabs/nf-lamin), as one artifact per dataset store.

Building a dataset takes two steps, each a Nextflow run of this repository launched through Seqera from the `build/main` branch:

1. A load script runs the loader of a source and publishes the raw datasets, e.g. `scrnaseq/cellxgene_census/hlca/raw.anndata.zarr`.
2. `process_datasets.py` processes the latest raw dataset of every id in a modality and publishes the common datasets, e.g. `scrnaseq/cellxgene_census/hlca/dataset.anndata.zarr`. It reads each raw dataset by its artifact uid, so LaminDB links the two.

| Script | Workflow | Datasets |
|---|---|---|
| `cellxgene_census.py` | `scrnaseq/load_cellxgene_census` | 8 CELLxGENE Census datasets (Census 2025-11-08) |
| `openproblems_op3.py` | `scrnaseq/load_openproblems_op3` | OP3 |
| `openproblems_neurips2021_bmmc.py` | `multimodal/load_openproblems_neurips2021_bmmc` | NeurIPS 2021 CITE-seq and Multiome |
| `openproblems_neurips2022_pbmc.py` | `multimodal/load_openproblems_neurips2022_pbmc` | NeurIPS 2022 CITE-seq and Multiome |
| `process_datasets.py` | `scrnaseq/process_scrnaseq`, `multimodal/process_multimodal` | the latest raw datasets of a modality |

The spatial scripts (`tenx_visium.sh`, `zenodo.sh`, `zenodo_slidetags.sh`) still use the old setup.

## Running the scripts

Run them from the root of the repository. Each script lists its dependencies, so `uv` installs them:

```bash
# load
uv run scripts/create_resources/cellxgene_census.py
uv run scripts/create_resources/openproblems_op3.py
uv run scripts/create_resources/openproblems_neurips2021_bmmc.py
uv run scripts/create_resources/openproblems_neurips2022_pbmc.py

# once the load runs have finished: process
uv run scripts/create_resources/process_datasets.py scrnaseq
uv run scripts/create_resources/process_datasets.py multimodal

# once the process runs have finished: annotate
uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --key-prefix scrnaseq/
uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --key-prefix multimodal/
```

Every script takes these options:

* `--ids a,b`: only launch these datasets, e.g. `--ids cellxgene_census/dkd`.
* `--dry-run`: write the Nextflow config and parameters and print where they are, without launching.
* `--instance`: the Lamin instance to publish to. It defaults to `$LAMIN_INSTANCE` or `openproblems/data-dev`; any other instance also needs `ALLOW_INSTANCE=1`.

Every dataset is a run of its own, since nf-lamin links each artifact a run publishes to every artifact the run read. Datasets marked `large` (Tabula Sapiens, the Human Retina Cell Atlas, OP3 and the NeurIPS 2022 Multiome) get up to 450 GB of memory and 48 hours per step ([`src/config/large_datasets.config`](../../src/config/large_datasets.config)).

## Setup

* `tw` is logged in to the Seqera workspace with the `denbi_bibigrid_cpu` compute environment.
* The Seqera workspace has a secret `LAMIN_BOT_API_KEY`, the API key of a Lamin user that can write to the instance and read [`openproblems/private`](https://lamin.ai/openproblems/private), which holds the NeurIPS 2022 source files.
* `build/main` is up to date with the commit you want to run.

The Nextflow config is [`common/nextflow_helpers/labels_denbi.config`](../../common/nextflow_helpers/labels_denbi.config), then [`src/config/lamin_denbi.config`](../../src/config/lamin_denbi.config), which sets up nf-lamin, makes a failed dataset fail the run, and keeps the private inputs out of the public instance.

## Versions

Publishing a dataset again with the same content gives the same artifact. Publishing it with different content creates a new version, but nf-lamin currently overwrites the files of the previous version ([laminlabs/nf-lamin#78](https://github.com/laminlabs/nf-lamin/issues/78)), so only rerun a dataset on purpose until that is resolved.
