import sys
import os
import pytest
import anndata as ad
import numpy as np

## VIASH START
meta = {
    "resources_dir": "resources_test/common"
}
## VIASH END

input_path = f"{meta['resources_dir']}/cxg_mouse_pancreas_atlas/dataset.anndata.zarr"
input = ad.read_zarr(input_path)

def test_even_sampling(run_component):
    output_path = "output.zarr"
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--even",
        "--seed", "123",
        "--n_obs", "100",
        "--n_vars", "120"
    ])

    # Checking whether file exists
    assert os.path.exists(output_path), "Output file not found"

    # Check that test output fits expected API
    output = ad.read_zarr(output_path)

    assert output.n_obs <= 100, "n_obs should be <= 100"
    assert output.n_vars <= 120, "n_vars should be <= 100"


def test_keep_functionality(run_component):
    output_path = "output.zarr"

    # take categories and well-expressed genes from the test data
    cell_types = input.obs["cell_type"].value_counts().index[:2].tolist()
    batches = input.obs["batch"].value_counts().index[:2].tolist()
    expression = np.asarray(input.layers["counts"].sum(axis=0)).ravel()
    keep_features = input.var_names[np.argsort(-expression)[:3]].tolist()

    run_component([
        "--input", input_path,
        "--keep_cell_type_categories", ";".join(cell_types),
        "--keep_batch_categories", ";".join(batches),
        "--keep_features", ";".join(keep_features),
        "--output", output_path,
        "--seed", "123"
    ])

    # Checking whether file exists
    assert os.path.exists(output_path), "Output file not found"

    # Check that test output fits expected API
    output = ad.read_zarr(output_path)

    assert output.n_obs <= 500, "n_obs should be <= 500"
    assert output.n_vars <= 500, "n_vars should be <= 500"
    for feat in keep_features:
        assert feat in output.var_names, f"{feat} should be in output.var_names"
    assert set(output.obs["cell_type"]) <= set(cell_types), "only the kept cell types should remain"
    assert set(output.obs["batch"]) <= set(batches), "only the kept batches should remain"

if __name__ == '__main__':
    sys.exit(pytest.main([__file__, "--capture=no"], plugins=["viashpy"]))
