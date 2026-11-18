import shutil
import subprocess
from pathlib import Path

import anndata as ad

## VIASH START
meta = {
  "executable": "target/executable/validation/validate_dataset/validate_dataset",
  "resources_dir": "resources_test/common",
  "temp_dir": "/tmp",
}
## VIASH END

resources = Path(meta["resources_dir"]) / "cxg_mouse_pancreas_atlas"
raw = resources / "raw.anndata.zarr"
dataset = resources / "dataset.anndata.zarr"


def validate(*args):
  return subprocess.run([meta["executable"], *args], capture_output=True, text=True)


print(">> Valid raw and processed datasets pass", flush=True)
for args in [
  ["--input_raw", raw],
  ["--input_dataset", dataset],
  ["--input_raw", raw, "--input_dataset", dataset],
]:
  out = validate(*map(str, args))
  assert out.returncode == 0, out.stdout + out.stderr

print(">> A raw dataset is not a processed dataset", flush=True)
out = validate("--input_dataset", str(raw))
assert out.returncode != 0
assert "is missing slot" in out.stdout + out.stderr, out.stdout + out.stderr

print(">> A raw dataset without a required obs column fails", flush=True)
adata = ad.read_zarr(raw)
del adata.obs["cell_type_ontology_term_id"]
broken = Path(meta["temp_dir"]) / "broken.anndata.zarr"
shutil.rmtree(broken, ignore_errors=True)
adata.write_zarr(broken)
out = validate("--input_raw", str(broken))
assert out.returncode != 0
assert "cell_type_ontology_term_id" in out.stdout + out.stderr, out.stdout + out.stderr

print(">> Passing no input fails", flush=True)
assert validate().returncode != 0

print(">> All tests passed", flush=True)
