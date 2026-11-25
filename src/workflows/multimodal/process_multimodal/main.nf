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

    // one AnnData per modality, so the single-modality components can process them
    | split_mudata.run(
      fromState: ["input": "input"],
      toState: [
        "split_rna": "output_rna",
        "split_prot": "output_prot",
        "split_atac": "output_atac"
      ]
    )

    // normalise each modality in the way that suits it
    | log_cp.run(
      fromState: ["input": "split_rna"],
      args: [normalization_id: "log_cp10k", n_cp: 10000],
      toState: ["normalized_rna": "output"]
    )
    | prot_clr.run(
      runIf: { id, state -> state.split_prot },
      fromState: ["input": "split_prot"],
      args: [normalization_id: "prot_clr"],
      toState: ["normalized_prot": "output"]
    )
    | atac_tfidf.run(
      runIf: { id, state -> state.split_atac },
      fromState: ["input": "split_atac"],
      args: [normalization_id: "atac_tfidf"],
      toState: ["normalized_atac": "output"]
    )

    // embedding and highly variable features per modality
    | svd.run(
      key: "svd_rna",
      fromState: ["input": "normalized_rna"],
      args: [output_mod2: null],
      toState: ["svd_rna": "output"]
    )
    | hvg.run(
      key: "hvg_rna",
      fromState: ["input": "svd_rna"],
      toState: ["hvg_rna": "output"]
    )
    | svd.run(
      key: "svd_prot",
      runIf: { id, state -> state.normalized_prot },
      fromState: ["input": "normalized_prot"],
      args: [output_mod2: null],
      toState: ["svd_prot": "output"]
    )
    | hvg.run(
      key: "hvg_prot",
      runIf: { id, state -> state.normalized_prot },
      fromState: ["input": "svd_prot"],
      toState: ["hvg_prot": "output"]
    )
    | svd.run(
      key: "svd_atac",
      runIf: { id, state -> state.normalized_atac },
      fromState: ["input": "normalized_atac"],
      args: [output_mod2: null],
      toState: ["svd_atac": "output"]
    )
    | hvg.run(
      key: "hvg_atac",
      runIf: { id, state -> state.normalized_atac },
      fromState: ["input": "svd_atac"],
      toState: ["hvg_atac": "output"]
    )

    | combine_mudata.run(
      fromState: [
        "input_raw": "input",
        "input_rna": "hvg_rna",
        "input_prot": "hvg_prot",
        "input_atac": "hvg_atac"
      ],
      toState: ["output_dataset": "output"]
    )

    // fail before publishing anything that does not match the common multimodal format
    | validate_dataset.run(
      fromState: ["input_multimodal_dataset": "output_dataset"],
      toState: { id, output, state -> state }
    )

    | setState(["output_dataset"])

  emit:
  output_ch
}
