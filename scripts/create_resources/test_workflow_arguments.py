"""Check the parameters of every dataset against the arguments of its workflow.

Viash refuses a `param_list` entry with a key that is not an argument of the workflow, so
this catches a bad parameter before a launch does. It reads the built workflow configs, so
build them first with `viash ns build --setup cachedbuild`.

    uv run --with pytest --with pyyaml pytest scripts/create_resources
"""

import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
import cellxgene_census  # noqa: E402
import helpers  # noqa: E402
import openproblems_neurips2021_bmmc  # noqa: E402
import openproblems_neurips2022_pbmc  # noqa: E402
import openproblems_op3  # noqa: E402
import process_datasets  # noqa: E402

UID = "abcdefghijklmnop0000"


def workflow_arguments(workflow):
  config = helpers.REPO_ROOT / "target/nextflow/workflows" / workflow / ".config.vsh.yaml"
  if not config.exists():
    pytest.skip(f"{config} is not built")
  config = yaml.safe_load(config.read_text())
  arguments = config.get("arguments", []) + [a for group in config.get("argument_groups", []) for a in group["arguments"]]
  return {argument["name"].removeprefix("--") for argument in arguments}


@pytest.mark.parametrize("workflow,module", [
  ("scrnaseq/load_cellxgene_census", cellxgene_census),
  ("scrnaseq/load_openproblems_op3", openproblems_op3),
  ("multimodal/load_openproblems_neurips2021_bmmc", openproblems_neurips2021_bmmc),
  ("multimodal/load_openproblems_neurips2022_pbmc", openproblems_neurips2022_pbmc),
])
def test_load_params_are_arguments(workflow, module):
  arguments = workflow_arguments(workflow) | {"id"}
  for dataset in helpers.split_by_size(module.DATASETS)[0] + helpers.split_by_size(module.DATASETS)[1]:
    assert set(dataset) <= arguments, f"{dataset['id']}: {set(dataset) - arguments}"
  assert set(module.SHARED) <= arguments, set(module.SHARED) - arguments


@pytest.mark.parametrize("workflow,key", [
  ("scrnaseq/process_scrnaseq", "scrnaseq/cellxgene_census/hlca/raw.anndata.zarr"),
  ("multimodal/process_multimodal", "multimodal/openproblems_neurips2021/bmmc_cite/raw.mudata.zarr"),
])
def test_process_params_are_arguments(workflow, key):
  arguments = workflow_arguments(workflow) | {"id"}
  param = process_datasets.raw_key_to_param(key, UID, "openproblems/data-dev")
  assert set(param) <= arguments, set(param) - arguments
