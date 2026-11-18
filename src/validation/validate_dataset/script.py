import yaml
from openproblems.project.component_tests.run_and_check_output import check_format

## VIASH START
par = {
  "input_raw": "resources_test/common/cxg_mouse_pancreas_atlas/raw.anndata.zarr",
  "input_dataset": "resources_test/common/cxg_mouse_pancreas_atlas/dataset.anndata.zarr",
}
meta = {
  "config": "target/executable/validation/validate_dataset/.config.vsh.yaml",
}
## VIASH END

# read the file formats from the built config, in which viash has resolved the merges
with open(meta["config"]) as f:
  config = yaml.safe_load(f)
arguments = [
  arg
  for group in config.get("argument_groups", [])
  for arg in group["arguments"]
]

# check every input that was passed
inputs = [arg for arg in arguments if par.get(arg["name"].lstrip("-")) is not None]
assert inputs, "Pass at least one of " + ", ".join(arg["name"] for arg in arguments)

for arg in inputs:
  clean_name = arg["name"].lstrip("-")
  print(f">> Checking {clean_name} against '{arg['label']}'", flush=True)
  check_format({**arg, "clean_name": clean_name, "value": par[clean_name]})

print(">> All inputs match their format", flush=True)
