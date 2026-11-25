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

    | openproblems_neurips2022_pbmc.run(
      fromState: [
        "input_mod1": "input_mod1",
        "input_mod2": "input_mod2",
        "mod1": "mod1",
        "mod2": "mod2",
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
      fromState: ["input_multimodal_raw": "output_raw"],
      toState: { id, output, state -> state }
    )

    | setState(["output_raw"])

  emit:
  output_ch
}
