"""Launch dataset workflows on the de.NBI cluster through Seqera, publishing the datasets they
create to LaminDB with nf-lamin.

The Nextflow config is common/nextflow_helpers/labels_denbi.config, then
src/config/lamin_denbi.config, then the Lamin instance, and for large datasets
src/config/large_datasets.config. nf-lamin records each run in LaminDB, so the scripts
themselves are not tracked.
"""

import argparse
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INSTANCE = "openproblems/data-dev"
PREFIXES = ("scrnaseq", "multimodal")


@dataclass
class Launch:
  config: Path
  params: Path
  command: list[str]


def resolve_instance(instance=None):
  """The instance to publish to; anything but openproblems/data-dev needs ALLOW_INSTANCE=1."""
  instance = instance or os.environ.get("LAMIN_INSTANCE") or DEFAULT_INSTANCE
  if instance != DEFAULT_INSTANCE and os.environ.get("ALLOW_INSTANCE") != "1":
    raise ValueError(f"Refusing to publish to '{instance}', set ALLOW_INSTANCE=1 if you mean it")
  return instance


def launch(workflow, params, prefix, *, large=False, labels=None, instance=None, dry_run=False, work_dir=None):
  """Launch `target/nextflow/workflows/<workflow>/main.nf` with `params` on de.NBI.

  The datasets are published to `lamin://<instance>?prefix=<prefix>`. With `dry_run`, the
  config and params are written and the command is returned, but nothing is launched.
  """
  if prefix not in PREFIXES:
    raise ValueError(f"prefix must be one of {PREFIXES}, got '{prefix}'")
  instance = resolve_instance(instance)
  work_dir = Path(work_dir or tempfile.mkdtemp())

  parts = [
    (REPO_ROOT / "common/nextflow_helpers/labels_denbi.config").read_text(),
    (REPO_ROOT / "src/config/lamin_denbi.config").read_text(),
    f"lamin {{ instance = '{instance}' }}\n",
  ]
  if large:
    parts.append((REPO_ROOT / "src/config/large_datasets.config").read_text())
  config = work_dir / "nextflow.config"
  config.write_text("\n".join(parts))

  params_file = work_dir / "params.yaml"
  params_file.write_text(yaml.safe_dump({**params, "publish_dir": f"lamin://{instance}?prefix={prefix}"}, sort_keys=False))

  command = [
    "tw", "launch", "https://github.com/openproblems-bio/datasets.git",
    "--revision", "build/main",
    "--pull-latest",
    "--main-script", f"target/nextflow/workflows/{workflow}/main.nf",
    "--workspace", "53907369739130",
    "--compute-env", "denbi_bibigrid_cpu",
    "--params-file", str(params_file),
    "--config", str(config),
    "--labels", labels or f"{workflow.split('/')[-1]},dataset_build",
    "--workspace-secrets", "LAMIN_BOT_API_KEY",
  ]

  ids = ", ".join(entry["id"] for entry in params.get("param_list", []))
  print(f"{workflow} -> lamin://{instance}?prefix={prefix}{' (large)' if large else ''}: {ids}")
  print(f"  config: {config}\n  params: {params_file}")
  if not dry_run:
    subprocess.run(command, check=True)
  return Launch(config=config, params=params_file, command=command)


def select(datasets, ids):
  """The datasets with the given ids, or all of them if `ids` is None."""
  if ids is None:
    return datasets
  known = {dataset["id"] for dataset in datasets}
  unknown = [id for id in ids if id not in known]
  if unknown:
    raise ValueError(f"Unknown dataset ids: {', '.join(unknown)}")
  return [dataset for dataset in datasets if dataset["id"] in ids]


def split_by_size(datasets):
  """Split datasets into regular and large ones (`large: True`), dropping the flag."""
  regular = [{k: v for k, v in d.items() if k != "large"} for d in datasets if not d.get("large")]
  large = [{k: v for k, v in d.items() if k != "large"} for d in datasets if d.get("large")]
  return regular, large


def parse_args(description):
  parser = argparse.ArgumentParser(description=description)
  parser.add_argument("--ids", type=lambda value: value.split(","), help="comma-separated dataset ids to launch (default: all)")
  parser.add_argument("--instance", help=f"Lamin instance to publish to (default: $LAMIN_INSTANCE or {DEFAULT_INSTANCE})")
  parser.add_argument("--dry-run", action="store_true", help="write the config and params, but don't launch")
  return parser.parse_args()


def launch_datasets(workflow, datasets, shared, prefix, args):
  """Launch the selected datasets, the large ones separately with more memory and time."""
  regular, large = split_by_size(select(datasets, args.ids))
  for group, is_large in [(regular, False), (large, True)]:
    if group:
      launch(workflow, {"param_list": group, **shared}, prefix, large=is_large, instance=args.instance, dry_run=args.dry_run)
