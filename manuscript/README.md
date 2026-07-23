# Manuscript analysis (phontrast paper)

Reproducibility code for the laboratory-phonology **manuscript** (distinct from
the frozen 2026 conference bundle on the `labphon_2026` branch). This directory
is tracked on `main` but excluded from the package build (`.Rbuildignore`), and
its `data/` and `outputs/` subfolders are git-ignored so raw research data and
generated tables stay local.

## `rerun_paper_metrics_phontrast210.R`

Regenerates the three metric tables with **phontrast 2.1.0**, reporting the
corrected Monte-Carlo estimator (`method = "mc"`, the 2.1.0 default) alongside
`method = "legacy"` (the estimator the conference poster tables used). Run from
the repository root:

```bash
Rscript manuscript/rerun_paper_metrics_phontrast210.R
```

Outputs land in `manuscript/outputs/`:

- `pb52_pairwise_phontrast210.csv` — PB52 F1/F2, 45 vowel pairs
- `mfcc13_pairwise_phontrast210.csv` — SBCSAE 13-MFCC, 351 vowel pairs

Each block is skipped with a message if its input is missing, so you can run
either on its own.

## Inputs (place under `manuscript/data/`)

| file | source |
| --- | --- |
| `rerun_mfcc_sampled_tokens_20260624.csv` | committed on the `labphon_2026` branch under `analysis/labphon_2026/data/` |
| PB52 F1/F2 | loaded from the `phonTools` package — no file needed (`install.packages("phonTools")`) |

## Estimator notes for the manuscript

- **Headline result is estimator-robust.** JSD has the strongest correspondence
  with KDE overlap under both estimators. The advantage is a high-dimensional
  effect: in 2-D PB52 F1/F2 every metric corresponds nearly perfectly (all
  Spearman ρ ≈ 0.98–1.00), whereas in SBCSAE 13-MFCC JSD pulls clearly ahead
  (ρ ≈ 0.95 mc / 1.00 legacy, vs Pillai ≈ 0.88 and Bhattacharyya ≈ 0.79).
- **`mc` is primary** (consistent estimate of the continuous JSD; 2.1.0 default)
  with `legacy` reported as a robustness column.
- The manuscript's vowel categories are large (PB52 E/I: 152 + 152 tokens; each
  SBCSAE vowel ≈ 500), where `mc` does not floor and its finite-sample null bias
  is negligible. (The small-n flooring that motivated the 2.1.0 fix appears only
  well below these sizes.)
