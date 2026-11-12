"""Shared steps of the test resource notebooks in this directory.

The notebooks describe the datasets and how they are processed; the mechanics of
running the workflows and storing the results in LaminDB live here.
"""

import json
import os
import subprocess
import tempfile
from pathlib import Path

import lamindb as ln

INSTANCE = "openproblems/data-dev"
PROJECT = "openproblems-bio/datasets"
REPO_URL = "https://github.com/openproblems-bio/datasets"


def git(*args):
  return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout.strip()


REPO_ROOT = Path(git("rev-parse", "--show-toplevel"))


def check_instance():
  """Make sure the .lamin dev-dir in the repo root selected the right instance."""
  instance = ln.setup.settings.instance.slug
  if instance != INSTANCE:
    raise RuntimeError(f"Connected to '{instance}' instead of '{INSTANCE}'")


def link_to_github(notebook):
  """Point the transform of the current run to the notebook on GitHub.

  Only done when the notebook is committed and pushed, so the link shows the code that ran.
  """
  commit = git("rev-parse", "HEAD")
  if git("status", "--porcelain", "--", str(REPO_ROOT / notebook)) or not git("branch", "-r", "--contains", commit):
    print(f"Not linking to GitHub: {notebook} has uncommitted or unpushed changes")
    return
  transform = ln.context.transform
  transform.reference = f"{REPO_URL}/blob/{commit}/{notebook}"
  transform.reference_type = "url"
  transform.save()


def run_workflow(workflow, params, highmem=None):
  """Run a workflow of this repository, return the directory with its outputs.

  `highmem` overrides the memory of processes labelled `highmem` (e.g. "64.GB"),
  for workflows that load large source files.
  """
  output_dir = Path(tempfile.mkdtemp())
  params_file = output_dir / "params.json"
  params_file.write_text(json.dumps({**params, "publish_dir": str(output_dir)}))
  configs = ["-c", "common/nextflow_helpers/labels_ci.config"]
  if highmem is not None:
    memory_config = output_dir / "memory.config"
    memory_config.write_text(f"process {{\n  withLabel: highmem {{ memory = {highmem} }}\n}}\n")
    configs += ["-c", str(memory_config)]
  subprocess.run(
    [
      "nextflow", "run", ".",
      "-main-script", f"target/nextflow/workflows/{workflow}/main.nf",
      "-profile", "docker",
      *configs,
      "-params-file", str(params_file),
    ],
    cwd=REPO_ROOT,
    env={**os.environ, "NXF_SYNTAX_PARSER": "v1"},
    check=True,
  )
  return output_dir


def save_collection(output_dir, files, key, description):
  """Store each file as an artifact under `key` and group them in a collection with that key.

  Earlier versions keep their files, so a pinned version of the collection stays downloadable.
  """
  artifacts = [
    ln.Artifact(
      output_dir / file,
      key=f"{key}/{file}",
      description=f"{description}: {file}",
      overwrite_versions=False,
    ).save()
    for file in files
  ]
  collection = ln.Collection(artifacts, key=key, description=description).save()
  print(f"lamin://{INSTANCE}/collection/{collection.uid}")
  return collection
