# Test fixtures

## Peterson and Barney (1952): the published calibration of `rank_contrasts()`

`test-pb52-calibration.R` checks that `rank_contrasts()` reproduces, on the 45
Peterson--Barney F1 x F2 vowel pairs, the ranking-protocol outcome reported in
Berry, "Estimand or estimator? Comparing vowel overlap measures against a known
ground truth" (JASA, under review, JASA-14548R1), Sec. VII A.

- `pb52_f1f2.csv`: `vowel,f1,f2`, 1,520 tokens (ten vowels, 152 tokens each;
  X-SAMPA labels `3'`, `A`, `E`, `I`, `O`, `U`, `V`, `i`, `u`, `{`). The
  de-identified slice the paper's replication code reads (no speaker, sex, or
  repetition columns). Copied byte for byte from
  `projects/p-phontrast-methods/analysis/matched_bc/output/pb52_f1f2_deident.csv`
  in `berrygrant/lab-agents` (sha256 `bea041eb...18628ab`).
- `pb52_paper_per_pair.csv`: the paper's per-pair values, one row per pair:
  step 1 (`step1_sqrt_jsd`, `step1_pillai`, `step1_kde_overlap`; plug-in
  bandwidth), the ceiling (`step5_saturated`, sqrt(JSD) >= 0.99), the
  percentile ranks and flag over the 22 pairs below the ceiling (`*_Bprime`),
  and the bandwidth check (`step4_sqrt_jsd_x0p5`, `step4_sqrt_jsd_x2`,
  `step4_rank_shift_Bprime`, `step4_set_aside_Bprime`). A column subset of
  `projects/p-phontrast-methods/analysis/sim/output/procedure/procedure_per_pair_PB52.csv`
  in `berrygrant/lab-agents` (sha256 `5da96729...6457d11`), with the IPA
  display columns dropped so the file is ASCII.

The F1/F2 values are published measurements (Peterson, G. E., & Barney, H. L.
(1952). Control methods used in a study of the vowels. *JASA*, 24(2), 175--184,
<https://doi.org/10.1121/1.1906875>) and circulate as a standard public data
set, for example as `phonTools::pb52`.
