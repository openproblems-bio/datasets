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

    | log_cp.run(
      fromState: ["input": "input"],
      args: [normalization_id: "log_cp10k", n_cp: 10000],
      toState: ["output_normalized": "output"]
    )

    | hvg.run(
      fromState: ["input": "output_normalized"],
      toState: ["output_hvg": "output"]
    )

    | pca.run(
      fromState: ["input": "output_hvg"],
      toState: ["output_pca": "output"]
    )

    | knn.run(
      fromState: ["input": "output_pca"],
      toState: ["output_dataset": "output"]
    )

    // fail before publishing anything that does not match the common dataset format
    | validate_dataset.run(
      fromState: ["input_dataset": "output_dataset"],
      toState: { id, output, state -> state }
    )

    | setState([
      "output_dataset",
      "output_normalized",
      "output_hvg",
      "output_pca"
    ])

  emit:
  output_ch
}
