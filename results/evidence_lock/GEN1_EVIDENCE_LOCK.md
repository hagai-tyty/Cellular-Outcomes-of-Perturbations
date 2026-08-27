# CellFate-Rx Generation 1 — EVIDENCE LOCK

```text
  GEN1_EVIDENCE_LOCKED

  lock digest   99c35793162aaa0e02f681cfaf4d9488492bb712a567e29a061d4886287489e0
  artifacts     54
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
  code          785b9811cee22fe6b28173dbbe7bd109d519c81cba2592aa4410c6ed7a4cc7b9  experiments/build_stage22_prospective_benchmarks.py
  code          f81dc9714d01db01c98ef51a7f6b351b7e4f261a083d1947c732abc103181542  experiments/run_gen1_evidence_lock.py
  code          a0ae156176e68561175453c875036193e9719e0e30131e8ecb2e702994908d44  experiments/run_stage23_2h_confirmation.py
  code          ce43d831b0c7585a226a61c340c95660fa14e67ab2ea0fa3714d822c64e8b7f7  experiments/run_stage23_learnability_gate.py
  code          1649c9d41fd8304bc41353433b8970b16c555c40b64b69b4526368f0e5e8363a  experiments/run_stage24_gen1_tool.py
  code          45ea02f87358a0471d1535bfe1b1a4c7c381eefbb3a3b5986165128301fccc45  experiments/run_stage25_ranking.py
  code          952fb4493b00a620a276bfb98a66140451e4fe76af1e1c56ec0af270b667fd8e  experiments/run_stage26_scope_lock.py
  protocol      357c873a9ce04ab089f22c3e29ebcf9b67163a3b149745f8153f5ae362de9396  plans/(newer)practical plans/GEN1_EVIDENCE_LOCK_V1.md
  records       e9122a9343af3060285d6ddb105ec072e2b18c3ba0d0303d9d4982284bb94667  plans/(newer)practical plans/RECORDs/stage_22_RECORD.md
  records       5f30b143f7ea1cc3652b0a9d9af35ef4682ed72d4f1adc3a6e2aed0a56990df7  plans/(newer)practical plans/RECORDs/stage_23_2H_RECORD.md
  records       860e255415154e728f8eb41dcb1e6e324ce0b9a8b3fb6cf4a4cab44668f05b43  plans/(newer)practical plans/RECORDs/stage_23_5_RECORD.md
  records       beb5918f5198bc528cdebdcc11bfc50c5fa419d477d87e9230312b74195e2b20  plans/(newer)practical plans/RECORDs/stage_23_RECORD.md
  records       5c7e213bdbac7217df8d26adbc72d309b1ffbf6f64fa5aa430f036c500a787a8  plans/(newer)practical plans/RECORDs/stage_24_POSTFREEZE_ADDENDUM.md
  records       de641a7f8150476733f6c245965f1d78d34fba48597a414ac092801c1c284e65  plans/(newer)practical plans/RECORDs/stage_24_RECORD.md
  records       9d71abbc8ac6795f7cb623c15ad536f02caa44821a80b9921b2b36ee3fceb109  plans/(newer)practical plans/RECORDs/stage_25_RECORD.md
  records       2a3fa981d97e0e44c9500c8374e9bc6c22e402ce7d0b0854d53f50355a53f1c6  plans/(newer)practical plans/RECORDs/stage_26_RECORD.md
  protocol      713c5509538dd5e704f54661c288f02f97dc8ba7a857422818747f5403251b47  plans/(newer)practical plans/STAGE_23_2_ROLE_A_CONFIRMATION_V5.md
  protocol      d7b874ac2801919a5ab9c1898a4ee2b882aa0ffc08b8a364be186e545cddfb9b  plans/(newer)practical plans/STAGE_23_2_V5_ADDENDUM_1_POWER_CURVE.md
  protocol      59f22e9ad2f8ccc40056e0a163ee5a926bb2a6f01cc7992d3f8e03e75b453e01  plans/(newer)practical plans/STAGE_23_5_GEN1_ROLE_B_SHIP_PLAN_V1.md
  protocol      f34a85055071e59be0d169cace92ba38892b49b5d9f22e4b19248b9b3535e58d  plans/(newer)practical plans/STAGE_26_KNOWN_TREATMENT_SCOPE_LOCK_V1.md
  benchmark     750e47623b9bd473ce9a642d551352883b84b25c274c584617be839e964eb08f  results/stage22_wm989_clones.csv
  supporting_role_A  f022d452c21af1e299ed9620ae16f498156dfb8b6119a64948e38412a255cf0a  results/stage23_2/stage23_2_handoff_to_stage24.json
  supporting_role_A  0983a3c5c204ede773408a1e03fe1ac75a2dd2874a811444fa3d3c40821efcf9  results/stage23_2h/stage23_2h_confirmation.json
  supporting_role_A  2c4a174ab1812cdc4ab8ddf8e76bc60bc5437e76095e81936d5cc82fecceef4d  results/stage23_2h/stage23_2h_power.json
  supporting_role_A  9e5f40d11f98cf0c5f1e6dc3dd11a411c1f1463c1c516d0bb19b4dfba54b6c4c  results/stage23_2h/stage23_2h_power_audit.json
  supporting_role_A  0ac07f05b066b7cfdae87dafe2ba8cee0f6a9aed4765eca5d7684d974d19f0f7  results/stage23_2h/stage23_2h_verdict.json
  protocol      720d939ce3ef278c1995ae30848175471e5ef7ba1dab06c0dbdd89e1a1c63f62  results/stage23_5_handoff_to_stage24.json
  protocol      7e063190e9f6af2bfe9156d70d4ae4d826c41d806fa2023a34e748064024556a  results/stage23_5_protocol.json
  benchmark     ba09c1933e3039f1b097b3a4dba6f2df967e4d31c1e3509b29e0a50fc31aa820  results/stage23_wm989_detection_oof.csv
  benchmark     e727c5782ac67ec6324219f4662dee72b00285dd0562435065941a9d831874aa  results/stage23_wm989_interaction_oof.csv
  predictions   abdc2999c0ee388fb2f2928318ff8a935921aaa72462d6d3c0ed64384e12e9a6  results/stage24/stage24_oof_for_stage25.csv
  tool          bda2c2de80187d0fd059da3c12be3f718422c191af4d941cc48748af94af59bc  results/stage24/stage24_w5_artifact.json
  tool          954cef7cff296d994d10bae3741ced1c8c9c5c7f97b13877a2916cf26b713315  results/stage24/stage24_w5_artifact.npz
  verdicts      ee8e2652a398ef91687104d3add5dc72c5e7c7cfb3adb7de0fb98d04c365d2e5  results/stage24/stage24f_tool_freeze.json
  tool          489f826266b8edd7cee17edec41998d0acb07d8c733ffe8e94fcb69bfc066464  results/stage24/tool/MODEL_CARD.md
  tool          a470d22bab529c3a6d3d0de1809cd48cbfbf9974f8f94adfe8dcc2eb9859f4e9  results/stage24/tool/example_clone_README.md
  tool          f4f203af44d28182ebd4fa4a3b2c2c53f7a6afb682486ddb7cfe2c3d28f35020  results/stage24/tool/example_clone_expression.npy
  tool          903db8d0e200ac615e8b27425730e2443029206259e20e5c91f8fd6cede3db8d  results/stage24/tool/example_clone_nuisance.txt
  tool          330cc7c8a68d77c0eabfef5da015f48c8969f0ce7cd1a4281a73af4028621c2e  results/stage24/tool/example_clones.csv
  tool          ceb798affe8d1a7e1578055eace0c293b47b1eedb5a9ed792800d5ca47a97aec  results/stage24/tool/io_schema.json
  verdicts      904992ca7790d10d13fdddc55d6f3602701f7a1ae0f122de7020caafbf329e62  results/stage24_handoff_to_stage25.json
  verdicts      c8c40f762203fb92f0375474647ddf5569c79a565fb8a1541bc1bc034b76b943  results/stage25/stage25_verdict.json
  verdicts      97fd48adfa07664769563da16094417f142087aeea985ad4c690dd7483f4b3ec  results/stage25/stage25a_observed.json
  limitations   90563ada20d3c9aec52f8c3425b3f3ed1181eaa0d27b39aed558248bcd31d4e3  results/stage26/GEN1_SCOPE_LIMIT.md
  verdicts      0a60efcb09dbc327a3ca690c432628b8365766f93ca32b56e84a328a42479ea7  results/stage26/stage26_verdict.json
  verdicts      90292bae5c0de7ab62b36b438156a7e66875ae5f66b20d82c7a1e3a21c84efd0  results/stage26_handoff_to_evidence_lock.json
  tool          4ec35e3862c592ff1c25bf390dd00ba3d0203777f0585963153491b048fabae3  src/cellfate/gen1_cli.py
  tool          f9194e0d5db6017dc6286539658abb1f2cb0c573fb1c52786b18a28451ed5c94  src/cellfate/gen1_predictor.py
  code          9575e8301a33d0017ed245816201c1a3016f2717a282d9495feed35a42a9ca2c  tests/test_gen1_evidence_lock.py
  code          fdc3939d028dae01440b748f485501707f03b59dd46fdffad2bf42984e3342ea  tests/test_gen1_predictor.py
  code          577ac5e73d7d4287d78bf993352da098c15a90fa0d71ff2423a557a4781c6bd1  tests/test_stage23_2h_confirmation.py
  code          7c1b9bb0f88eefa1075d3bb546a841b8431d4109ccc77143b4314676bcb3241c  tests/test_stage24_gen1_tool.py
  code          a880c0163b00ad60e5848c848e05501e10a6fc776207eab3890afefdca6f0d4d  tests/test_stage25_ranking.py
  code          66ea8bed1c42d611dfe42b1777004e3effc75e7c150f3c757142767cd7a98d16  tests/test_stage26_scope_lock.py
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
