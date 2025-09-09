workflow run_wf {
  take:
  input_ch

  main:
  output_ch = input_ch

    // extract the dataset metadata
    | extract_uns_metadata.run(
      fromState: [input: "input"],
      toState: [output: "output"]
    )

    | setState([
      "output",
    ])

  emit:
  output_ch
}
