# LabPhon 2026 Poster Reproducibility

This folder is the reproducible analysis bundle for the final LabPhon 2026
poster:

`poster/LabPhon2026_Poster_Berry_Final.pdf`

The poster PDF is A0 portrait. Its SHA-256 hash is:

`2fcd75b61204cdea5f9868a6b737951cade3fdb3f4fc6816c2c377e92b12c558`

## What This Folder Reproduces

Running the Python script below regenerates the poster-facing analysis tables
and figures from committed processed inputs:

```bash
python3 -m pip install -r analysis/labphon_2026/requirements-python.txt
python3 analysis/labphon_2026/reproduce_poster_figures.py
```

The script writes regenerated outputs to `analysis/labphon_2026/outputs/`,
including:

- `poster_main_metric_comparison_20260624.{png,pdf,svg}`
- `poster_overlap_estimator_audit_no_mmo_20260625.{png,pdf,svg}`
- `stanley_jsd_sample_size_robustness.{png,pdf,svg}`
- `final_consonant_jsd_generalization.{png,pdf,svg}`
- `poster_pb52_jsd_mmo_comparison_20260625.{png,pdf,svg}`
- `poster_pb52_mmo_20260625.{png,pdf,svg}`
- `audited_overlap_estimator_values_no_mmo_20260625.csv`
- `audited_overlap_estimator_metric_performance_no_mmo_20260625.csv`
- `audited_overlap_estimator_method_notes_no_mmo_20260625.csv`
- `poster_result_summary_20260625.csv`

## Inputs

Processed inputs are committed under `data/`.

- `rerun_pb52_pairwise_20260624.csv`: PB52 F1/F2 pairwise metric table.
- `rerun_mfcc_sampled_tokens_20260624.csv`: SBCSAE sampled 13-MFCC vowel tokens.
- `rerun_pairwise_metrics_mfcc13_20260624.csv`: SBCSAE 13-MFCC pairwise metric table.
- `stanley_jsd_sample_size_summary.csv`: sample-size simulation summary under complete overlap.
- `final_consonant_jsd_generalization.csv`: exploratory stop-voicing and fricative-place JSD summary.
- `pb52_E_I_mmo_*`: faithful PB52 MMO talking-point outputs and diagnostics.

The bundled data are processed analysis inputs, not raw audio/TextGrid files.
The poster findings can be reproduced from these committed inputs. Raw SBCSAE
audio/alignment extraction is intentionally not required for checking the final
poster numbers.

## Main Poster Claims Encoded Here

- JSD has the strongest correspondence with KDE-estimated vowel overlap in both
  PB52 F1/F2 space and SBCSAE 13-MFCC space.
- The KDE comparison is intentionally described as KDE-estimated overlap,
  because phonJSD also estimates continuous distributions using KDE.
- The overlap-estimator robustness graph does not include MMO. It compares KDE,
  kNN density, MVN Monte Carlo, PCA-2D grid, and PCA-2D convex-hull support.
- PCA-2D convex hull is a geometric support-overlap diagnostic, not a true
  distribution comparison.
- The MMO scripts are retained as a separate PB52 talking-point analysis.

## Optional Faithful PB52 MMO Rerun

The faithful PB52 MMO analysis is computationally heavier and requires `brms`
and a working Stan backend. The committed outputs are:

- `data/pb52_E_I_mmo_ba_summary.csv`
- `data/pb52_E_I_mmo_diagnostics.csv`
- `data/pb52_E_I_mmo_ba_draws.csv`
- `data/pb52_E_I_mmo_model_data.csv`

To rerun:

```bash
Rscript analysis/labphon_2026/run_pb52_mmo_brms.R --fit
Rscript analysis/labphon_2026/compare_pb52_jsd_mmo.R
python3 analysis/labphon_2026/reproduce_poster_figures.py
```

The PB52 MMO model used here is faithful to the available PB52 data:
speaker-normalized F1/F2, vowel + talker type + repetition, and speaker random
effects. PB52 does not contain lexical frequency or phonological-context
covariates, so those controls from the full MMO framework cannot be included.

## Cleanup Decision

The old `labphon_2026` branch contained exploratory analysis scripts and draft
outputs. Those were not merged into `main` before cleanup because `main` already
contains the current package release and tests, while the old scripts were
poster-specific and superseded by this final bundle. They remain recoverable in
the earlier branch history.
