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

    // a moving release would give a different dataset for the same parameters
    | map{ id, state ->
      if (!state.input_uri && !(state.census_version ==~ /\d{4}-\d{2}-\d{2}/)) {
        error "Dataset '${id}': --census_version must be the date of a Census release, got '${state.census_version}'"
      }
      [id, state]
    }

    // fetch data from the CELLxGENE Census
    | cellxgene_census.run(
      fromState: [
        "input_uri": "input_uri",
        "census_version": "census_version",
        "species": "species",
        "obs_value_filter": "obs_value_filter",
        "cell_filter_grouping": "cell_filter_grouping",
        "cell_filter_minimum_count": "cell_filter_minimum_count",
        "obs_batch": "obs_batch",
        "obs_batch_separator": "obs_batch_separator",
        "dataset_id": "id",
        "dataset_name": "dataset_name",
        "dataset_url": "dataset_url",
        "dataset_reference": "dataset_reference",
        "dataset_summary": "dataset_summary",
        "dataset_description": "dataset_description",
        "dataset_organism": "dataset_organism",
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
