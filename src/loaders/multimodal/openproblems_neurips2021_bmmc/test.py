from os import path
import subprocess
import mudata as md

input = "neurips2021_bmmc_cite.zarr"
mod1 = "GEX"
mod2 = "ADT"

output_file = "output.mudata.zarr"

input_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122%5Fopenproblems%5Fneurips2021%5Fcite%5FBMMC%5Fprocessed%2Eh5ad%2Egz"

# download input
print(">> Downloading input", flush=True)
out = subprocess.run(
  [
    "wget",
    "-O", input + ".gz",
    input_url,
  ],
  stderr=subprocess.STDOUT
)
# unzip input
print(">> Unzipping input", flush=True)
out = subprocess.run(
  [
    "gunzip",
    input + ".gz",
  ],
  stderr=subprocess.STDOUT
)

print(">> Running script", flush=True)
out = subprocess.run(
  [
    meta["executable"],
    "--input", input,
    "--mod1", mod1,
    "--mod2", mod2,
    "--output", output_file,
    "--dataset_id", "openproblems/neurips2021_bmmc",
    "--dataset_name", "BMMC (CITE-seq)",
    "--dataset_url", "http://foo.org",
    "--dataset_reference", "foo2000bar",
    "--dataset_summary", "A short summary.",
    "--dataset_description", "A couple of paragraphs worth of text.",
  ],
  stderr=subprocess.STDOUT
)

if out.stdout:
    print(out.stdout, flush=True)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.", flush=True)
    exit(out.returncode)

print(">> Checking whether the output exists", flush=True)
assert path.exists(output_file), "Output file does not exist"

print(">> Read output mudata", flush=True)
mdata = md.read_zarr(output_file)
print(mdata, flush=True)

print(">> Check that the output fits the expected API", flush=True)
assert set(mdata.mod) == {"rna", "prot"}, f"Expected modalities rna and prot, got {list(mdata.mod)}"
for name, mod in mdata.mod.items():
  assert mod.X is None, f".X is not empty in modality {name}"
  assert "counts" in mod.layers, f"'counts' not found in the layers of modality {name}"
  assert "cell_type" not in mod.obs.columns, f"cell_type is duplicated in the obs of modality {name}"
assert "cell_type" in mdata.obs.columns, "cell_type column not found in global obs"
assert "batch" in mdata.obs.columns, "batch column not found in global obs"
assert mdata.obs["organism_ontology_term_id"].unique().tolist() == ["NCBITaxon:9606"], "organism_ontology_term_id is not NCBITaxon:9606"
assert all(":" not in col for col in mdata.obs.columns), "global obs contains prefixed columns"
assert mdata.uns["dataset_id"] == "openproblems/neurips2021_bmmc"
assert mdata.uns["dataset_name"] == "BMMC (CITE-seq)"
assert mdata.uns["dataset_url"] == "http://foo.org"
assert list(mdata.uns["dataset_reference"]) == ["foo2000bar"]
assert mdata.uns["dataset_summary"] == "A short summary."
assert mdata.uns["dataset_description"] == "A couple of paragraphs worth of text."
assert "dataset_organism" not in mdata.uns

print(">> All tests passed successfully", flush=True)
