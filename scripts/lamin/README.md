# LaminDB scripts

`annotate_artifacts.py` annotates dataset artifacts in LaminDB with everything the OpenProblems website shows about a dataset: a schema with every slot of its file format in `src/api`, the ontology terms in `obs`, the `dataset_*` metadata, the structure of the store, cell counts per category, processing parameters, the modality and references with their DOI.

Run it after anything that published datasets, be it a Nextflow run or a script. Its dependencies are listed in the script, so `uv` installs them:

```bash
# everything a run created
uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --run <run uid>

# the latest version of every artifact under a key prefix
uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --key-prefix scrnaseq/

# specific artifacts, with a file format other than the one derived from the key
uv run scripts/lamin/annotate_artifacts.py --instance openproblems/data-dev --artifact <uid> --spec src/api/file_raw.yaml
```

The file format of an artifact follows from the end of its key (`raw.anndata.zarr` is a raw dataset, `dataset.anndata.zarr` a common dataset, `raw.mudata.zarr` and `dataset.mudata.zarr` their multimodal counterparts) and the modality from its start (`scrnaseq/...`); `--spec` and `--modality` override them. Add `--dry-run` to see what would be annotated.

For AnnData stores, only metadata, `obs`, `var` and `uns` are read, straight from where an artifact is stored, so this works for datasets of any size. MuData stores are annotated per modality, but lamindb 2.12 loads them into memory to validate them, so a large multimodal dataset needs the memory to hold it. Artifacts that already have the schema of their file format are skipped, unless you pass `--force`.
