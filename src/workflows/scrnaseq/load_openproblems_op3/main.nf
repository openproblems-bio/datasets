workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // fetch data from OP3 dataset
    | openproblems_op3.run(
      fromState: [
        "input": "input",
        "data_type": "data_type",
        "donor_id": "donor_id",
        "cell_type": "cell_type",
        "perturbation": "perturbation",
        "dataset_id": "id",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description"
      ],
      toState: ["output_raw": "output"]
    )

    // subsample if so desired
    | subsample.run(
      runIf: { id, state -> state.do_subsample },
      fromState: [
        "input": "output_raw",
        "n_obs": "n_obs",
        "n_vars": "n_vars",
        "keep_features": "keep_features",
        "keep_cell_type_categories": "keep_cell_type_categories",
        "keep_batch_categories": "keep_batch_categories",
        "even": "even",
        "seed": "seed"
      ],
      args: [output_mod2: null],
      toState: ["output_raw": "output"]
    )

    // fail before publishing anything that does not match the raw format
    | validate_dataset.run(
      fromState: ["input_raw": "output_raw"],
      toState: { id, output, state -> state }
    )

    | setState(["output_raw"])

  emit:
  output_ch
}
