"""Tests for helpers.py, without launching anything.

    uv run --with pytest --with pyyaml pytest scripts/create_resources
"""

import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).parent))
import helpers  # noqa: E402

DATASETS = [
  {"id": "a/small"},
  {"id": "a/big", "large": True},
  {"id": "b/small"},
]


def launch(tmp_path, **kwargs):
  args = {"workflow": "scrnaseq/load_cellxgene_census", "params": {"param_list": [{"id": "foo"}]}, "prefix": "scrnaseq"}
  return helpers.launch(**{**args, **kwargs}, dry_run=True, work_dir=tmp_path)


def test_config_order(tmp_path):
  config = launch(tmp_path).config.read_text()
  denbi = config.index("executor = 'slurm'")
  lamin = config.index("nf-lamin@")
  instance = config.index("instance = 'openproblems/data-dev'")
  assert denbi < lamin < instance
  assert "450.GB" not in config


def test_error_strategy(tmp_path):
  # retry when a task ran out of memory or its node went away (no or an odd exit code), and
  # let the other datasets finish when one fails, without ever ignoring a failure
  config = (helpers.REPO_ROOT / "src/config/lamin_denbi.config").read_text()
  strategy = next(line for line in config.splitlines() if line.strip().startswith("errorStrategy"))
  for condition in ["task.exitStatus in [137, 139]", "task.exitStatus == null", "task.exitStatus < 0", "task.exitStatus > 2100000000"]:
    assert condition in strategy
  assert "'finish'" in strategy
  assert "'ignore'" not in strategy and "'terminate'" not in strategy


def test_nf_lamin_pinned():
  # runs that start together install a new plugin version into the same shared folder on
  # /vol/scratch at once, and leave it broken
  config = (helpers.REPO_ROOT / "src/config/lamin_denbi.config").read_text()
  assert "id 'nf-lamin@1.0.1'" in config


def test_state_files_not_published():
  # Viash publishes a state yaml next to the outputs, which would end up in the instance's
  # storage without being registered
  config = (helpers.REPO_ROOT / "src/config/lamin_denbi.config").read_text()
  assert "withName: 'publishStatesProc' { publishDir = [enabled: false] }" in config


def test_large_config(tmp_path):
  config = launch(tmp_path, large=True).config.read_text()
  assert config.index("instance = 'openproblems/data-dev'") < config.index("450.GB")


def test_params_publish_to_lamin(tmp_path):
  params = yaml.safe_load(launch(tmp_path, params={"param_list": [{"id": "foo"}], "publish_dir": "s3://elsewhere"}).params.read_text())
  assert params["publish_dir"] == "lamin://openproblems/data-dev?prefix=scrnaseq"
  assert params["param_list"] == [{"id": "foo"}]


def test_tw_command(tmp_path):
  result = launch(tmp_path)
  command = result.command
  assert command[:3] == ["tw", "launch", "https://github.com/openproblems-bio/datasets.git"]
  assert command[command.index("--main-script") + 1] == "target/nextflow/workflows/scrnaseq/load_cellxgene_census/main.nf"
  assert command[command.index("--compute-env") + 1] == "denbi_bibigrid_cpu"
  assert command[command.index("--revision") + 1] == "build/main"
  assert command[command.index("--workspace-secrets") + 1] == "LAMIN_BOT_API_KEY"
  assert command[command.index("--params-file") + 1] == str(result.params)
  assert command[command.index("--config") + 1] == str(result.config)


def test_other_instance_refused(tmp_path, monkeypatch):
  monkeypatch.delenv("ALLOW_INSTANCE", raising=False)
  with pytest.raises(ValueError, match="ALLOW_INSTANCE"):
    launch(tmp_path, instance="openproblems/data")


def test_other_instance_allowed(tmp_path, monkeypatch):
  monkeypatch.setenv("ALLOW_INSTANCE", "1")
  params = yaml.safe_load(launch(tmp_path, instance="openproblems/data", prefix="multimodal").params.read_text())
  assert params["publish_dir"] == "lamin://openproblems/data?prefix=multimodal"


def test_unknown_prefix_refused(tmp_path):
  with pytest.raises(ValueError, match="prefix"):
    launch(tmp_path, prefix="spatial")


def test_select_by_id():
  assert helpers.select(DATASETS, None) == DATASETS
  assert helpers.select(DATASETS, ["a/big"]) == [{"id": "a/big", "large": True}]
  with pytest.raises(ValueError, match="c/missing"):
    helpers.select(DATASETS, ["a/big", "c/missing"])


def test_split_by_size():
  regular, large = helpers.split_by_size(DATASETS)
  assert regular == [{"id": "a/small"}, {"id": "b/small"}]
  assert large == [{"id": "a/big"}]
