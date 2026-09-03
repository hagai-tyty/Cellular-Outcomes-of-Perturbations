# CellFate-Rx Generation 1 — EVIDENCE LOCK

```text
  GEN1_EVIDENCE_LOCKED

  lock digest   9245e605f6272aa809858d3f32dbe55ed53864df90e0b55750ab0a7d577da400
  artifacts     62
  ship plan     8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48
```

The lock digest is a SHA-256 over `path  sha256` lines, sorted by path and LF-joined. It is one
number that names the entire Generation-1 evidence base, and it belongs in the manuscript.

**Two digest forms appear below and they are not in conflict.** The manifest hashes raw bytes,
which is what a file manifest must do. The frozen protocol identity `8da16fca...` is a
*canonical-LF* digest of the same ship plan, computed with CRLF normalised to LF so the protocol
has one identity on every platform. The ship plan therefore carries `59f22e9a...` as a file and
`8da16fca...` as a protocol, and both are checked.

## What is locked

```text
  protocol      1af392067f1d8283fbeb7211bbdfffe7dccda62fac811b31695ca18b8b30191e  CITATION.cff
  protocol      99a2884f6210ed76e569d2fc37bdbfcd0b0402b7fe972c2043fb2fb30a604c5a  environment_lock.txt
  code          dfd047eff9fb82d0303c24743416d5a3ddc524076976074adceb435cb5d704a7  experiments/build_stage22_prospective_benchmarks.py
  code          ff0bded1828ca59e06cdbad88dada45fa32b40cdeb8bcc583e4e1b275fb753fe  experiments/export_gen1_source_data.py
  code          75321aebe9671002620f2f635999b03754a9dd8fe6ec7511397de98faf32c2f6  experiments/make_gen1_figures.py
  code          87ee28a4a22256fda9103a8ac4fda8133073cc49a882a7c0b78dce017ae891b8  experiments/make_release_bundle.py
  code          cdbcfcbae793f4ed3de4083ddbdadf8c45e5100bf5e90438f73d315f5aa77dd4  experiments/run_gen1_evidence_lock.py
  code          a0801d4d2a91259b53e67c4705da0ccde7c822b73b9670cec9be1919614a0085  experiments/run_stage23_2h_confirmation.py
  code          ce43d831b0c7585a226a61c340c95660fa14e67ab2ea0fa3714d822c64e8b7f7  experiments/run_stage23_learnability_gate.py
  code          f04743532ef66cca47eb461cc4fdf304b53e0fbf113567ae3e7218f13b7b4ecc  experiments/run_stage24_gen1_tool.py
  code          2c06ff784c535647448b01545428d90aa7fccea5f5b8fda40567634dc386fd22  experiments/run_stage25_ranking.py
  code          1f207d5f4253a2aa6afc8aac349b6bb425922fa3938e72fd7046e3a0f6d9c511  experiments/run_stage26_scope_lock.py
  protocol      661cc2bf927cf070d1fb78f3bb936690303634b1d8de555a9a3bf29caa760390  plans/(newer)practical plans/GEN1_EVIDENCE_LOCK_V1.md
  records       e42cd1a4c132a3ae94e77a32bf67edfc3a0f69199ad0d7d776b10c1951f6c2dc  plans/(newer)practical plans/RECORDs/stage_22_RECORD.md
  records       bff1365349db9dba981f75e385cac54f4e446acfe9339527701425c2e8f057df  plans/(newer)practical plans/RECORDs/stage_23_2H_RECORD.md
  records       e42e57107fcd7abf5ec8dfa9c8dc209cfed637fbc99934089ffc6f99eded977f  plans/(newer)practical plans/RECORDs/stage_23_5_RECORD.md
  records       beb5918f5198bc528cdebdcc11bfc50c5fa419d477d87e9230312b74195e2b20  plans/(newer)practical plans/RECORDs/stage_23_RECORD.md
  records       5c7e213bdbac7217df8d26adbc72d309b1ffbf6f64fa5aa430f036c500a787a8  plans/(newer)practical plans/RECORDs/stage_24_POSTFREEZE_ADDENDUM.md
  records       8c3e9551c54b384e72a144c4fb1c99885dd2aa8eb8490687efd7a904b3ff0cac  plans/(newer)practical plans/RECORDs/stage_24_RECORD.md
  records       9d71abbc8ac6795f7cb623c15ad536f02caa44821a80b9921b2b36ee3fceb109  plans/(newer)practical plans/RECORDs/stage_25_RECORD.md
  records       2a3fa981d97e0e44c9500c8374e9bc6c22e402ce7d0b0854d53f50355a53f1c6  plans/(newer)practical plans/RECORDs/stage_26_RECORD.md
  protocol      b576c95aa1036470db1d041e2540f1ca66f4aee1e60e0f088edd34d1b4ba04fb  plans/(newer)practical plans/STAGE_23_2_ROLE_A_CONFIRMATION_V5.md
  protocol      f318623a14c203a0f7b0733d4c398ceb52296fad2231ec5d07773d757499f9f8  plans/(newer)practical plans/STAGE_23_2_V5_ADDENDUM_1_POWER_CURVE.md
  protocol      8da16fca0f84b5664f4668f86ed21530242be89020059d1c7ba98f22d7bced48  plans/(newer)practical plans/STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md
  protocol      f34a85055071e59be0d169cace92ba38892b49b5d9f22e4b19248b9b3535e58d  plans/(newer)practical plans/STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md
  benchmark     750e47623b9bd473ce9a642d551352883b84b25c274c584617be839e964eb08f  results/stage22_wm989_clones.csv
  supporting_role_A  f022d452c21af1e299ed9620ae16f498156dfb8b6119a64948e38412a255cf0a  results/stage23_2/stage23_2_handoff_to_stage24.json
  supporting_role_A  0983a3c5c204ede773408a1e03fe1ac75a2dd2874a811444fa3d3c40821efcf9  results/stage23_2h/stage23_2h_confirmation.json
  supporting_role_A  2c4a174ab1812cdc4ab8ddf8e76bc60bc5437e76095e81936d5cc82fecceef4d  results/stage23_2h/stage23_2h_power.json
  supporting_role_A  9e5f40d11f98cf0c5f1e6dc3dd11a411c1f1463c1c516d0bb19b4dfba54b6c4c  results/stage23_2h/stage23_2h_power_audit.json
  supporting_role_A  0ac07f05b066b7cfdae87dafe2ba8cee0f6a9aed4765eca5d7684d974d19f0f7  results/stage23_2h/stage23_2h_verdict.json
  protocol      f9fff3e6dfd957ba04387598ea8293b4c2777982aec86a97a39300d50daf084d  results/stage23_5_handoff_to_stage24.json
  protocol      c5feed957bfa8e31ee487ca1c44fcabb1286b615a6c0496997a26dd3a0c71454  results/stage23_5_protocol.json
  benchmark     ba09c1933e3039f1b097b3a4dba6f2df967e4d31c1e3509b29e0a50fc31aa820  results/stage23_wm989_detection_oof.csv
  benchmark     e727c5782ac67ec6324219f4662dee72b00285dd0562435065941a9d831874aa  results/stage23_wm989_interaction_oof.csv
  predictions   abdc2999c0ee388fb2f2928318ff8a935921aaa72462d6d3c0ed64384e12e9a6  results/stage24/stage24_oof_for_stage25.csv
  tool          2ecbc3c56bae2e9c9167068353f7ba4001646f74044be945c65823d6e36121ef  results/stage24/stage24_w5_artifact.json
  tool          954cef7cff296d994d10bae3741ced1c8c9c5c7f97b13877a2916cf26b713315  results/stage24/stage24_w5_artifact.npz
  verdicts      66f436e10519015ddff5ace9c3f91d8574f6f6e4a02a8514fdf8c979348bbdf3  results/stage24/stage24f_tool_freeze.json
  tool          f5d9b8334a1ec12e62dae5fb4d8cfd87dc8a177d365c2fbcfdf1c37f2593fd53  results/stage24/tool/MODEL_CARD.md
  tool          ae582114c8979180eef56b1570f61ef485b56b9e0a54b1b61b8c8f1687b7eeba  results/stage24/tool/example_clone_README.md
  tool          f4f203af44d28182ebd4fa4a3b2c2c53f7a6afb682486ddb7cfe2c3d28f35020  results/stage24/tool/example_clone_expression.npy
  tool          903db8d0e200ac615e8b27425730e2443029206259e20e5c91f8fd6cede3db8d  results/stage24/tool/example_clone_nuisance.txt
  tool          330cc7c8a68d77c0eabfef5da015f48c8969f0ce7cd1a4281a73af4028621c2e  results/stage24/tool/example_clones.csv
  tool          4e158b19f3e450628f661508992ef7d1ba7d83098b62dd60af3b9175e9edf021  results/stage24/tool/io_schema.json
  verdicts      e9910bd4c840dd2bed6247cc947781793cf9bed8a2179f0102d5675205355b13  results/stage24_handoff_to_stage25.json
  predictions   49ee6fb3bbe3a551d8164e622c1dab527813e3e0bf9ba51dff36ae8067ab2128  results/stage25/stage25_bootstrap_replicates.csv
  predictions   b0b4745727afb8a1ec9d3ab9f438592b3d514f08f7954958b18f51bbeee5dd5d  results/stage25/stage25_null_draws.csv
  verdicts      9355be0f3e94e82292f82f678b86bf624dd8fb11e9b6baecd37347c78f9e248d  results/stage25/stage25_verdict.json
  verdicts      00af0c8a33a061d609d808c0d02e4d710ccfa7baebe18f7c2edcb2d2757151f3  results/stage25/stage25a_observed.json
  limitations   e24081a49b13bc87e35ae95649d30d1717dc79e250c0b3c5110170ddc71ab6ea  results/stage26/GEN1_SCOPE_LIMIT.md
  verdicts      8e9888f01ae79884370aa55d96536ba9d62685db9cc457b0ad948e597de25295  results/stage26/stage26_verdict.json
  verdicts      48aa2a8aa4ffad2fcd7965cddac2450e0c4b8f3222873e7ea2affd1a7f4a2bb5  results/stage26_handoff_to_evidence_lock.json
  tool          4ec35e3862c592ff1c25bf390dd00ba3d0203777f0585963153491b048fabae3  src/cellfate/gen1_cli.py
  tool          6f7c84a2ab8beb14d64728f0241dcaee1bc39ba0abb1f879e52aca2410cf2944  src/cellfate/gen1_predictor.py
  code          1ca2450d1d58de9f0f7d5e60b8e7a4213df812aeebe25d272799d7c2461abe8b  tests/test_ci_portability.py
  code          004a48be15d1766a5e4bbb182c3ea3093e2ef069f614660d2a960c082c3d032b  tests/test_gen1_evidence_lock.py
  code          7b9a70e8b979e066b391548ffa59560d1e805e6866874c0f6e95114d8076e43c  tests/test_gen1_predictor.py
  code          7209b237c011ef8fa1498a5f046d204fbfa2600d63ce2cfefe326b5acf871223  tests/test_stage23_2h_confirmation.py
  code          0b71060fab5cc6ea04d935bfb72ce3cf3197f34e947439ed52fd7f405759c7c1  tests/test_stage24_gen1_tool.py
  code          a880c0163b00ad60e5848c848e05501e10a6fc776207eab3890afefdca6f0d4d  tests/test_stage25_ranking.py
  code          1f9951ae3def121b3c76465166e3678097fa0eb69829034bc80191b0c2ffe514  tests/test_stage26_scope_lock.py
```

## Chain of custody

An artifact's hash must be the same number everywhere a stage recorded it. If Stage 24 handed
Stage 25 a table hashed X and this lock hashes Y, the analysis did not consume what is being
locked.

```text
  out-of-fold table   identical across 4 independent records: True
  model artifact      identical across 4 independent records: True
  ship-plan digest    identical across 6 independent records: True
```

## The numbers this locks

```text
  eligible clones        892
  R(W1) / R(W4) / R(W5)  0.692654 / 0.692176 / 0.743781
  delta_RANK             +0.051605   CI95 [+0.037197, +0.065571]
  null p95               0.008672
  permutation            0 of 1000 draws reached the observed value
  p_perm                 0.000999   report as p < 0.001 (0 of 1,000), never as a point estimate
  delta_TOP1             +0.115471
  adversarial refusals   56 of 56
  design columns         309
```

## What this lock does NOT contain

```text
  stage24_w5_artifact.npz   44 MB, gitignored. A fresh clone does NOT contain it.
                            Its hash is locked. Rebuild:
                              python experiments/run_stage24_gen1_tool.py --stage 24c
  raw sequencing data       GSE279162 (WM989), GSE227151 (Rewind). Accessions are
                            locked; bytes are not vendored.
```

Naming a gap is not closing it. Both stay open.

## Verifying this lock

```text
  python experiments/run_gen1_evidence_lock.py --verify
```

It re-hashes every artifact and refuses if one has moved. Its ability to refuse is itself tested:
a one-bit flip and a deleted file must both be caught, on copies, before any lock is issued.

## What locking does not do

It grants no claim. It fixes what the existing claims are made of. No lock outcome reopens an
earlier stage, changes a recorded number, or authorizes new data, a new condition or a new model.
