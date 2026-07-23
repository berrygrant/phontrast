# Manuscript analysis (phontrast paper)

Reproducibility code for the laboratory-phonology **manuscript** — which pitches
`phontrast` as an R package for evaluating acoustic contrast and demonstrates its
utility for comparative analyses. Distinct from the frozen 2026 conference bundle
on the `labphon_2026` branch. This directory is tracked on `main` but excluded
from the package build (`.Rbuildignore`); its `data/` and `outputs/` subfolders
are git-ignored so raw research data and generated tables stay local.

## Scripts (run from the repository root, in order)

```bash
Rscript manuscript/rerun_paper_metrics_phontrast210.R   # 1. metric tables
Rscript manuscript/estimator_influence_audit.R          # 2. fairness audit table
Rscript manuscript/make_figures.R                       # 3. figures
```

1. **`rerun_paper_metrics_phontrast210.R`** — regenerates the metric tables with
   **phontrast 2.1.0**, reporting the corrected Monte-Carlo estimator
   (`method = "mc"`, the 2.1.0 default) alongside `method = "legacy"` (the
   estimator the conference poster used). Writes `pb52_pairwise_phontrast210.csv`
   (PB52 F1/F2, 45 pairs) and `mfcc13_pairwise_phontrast210.csv` (SBCSAE 13-MFCC,
   351 pairs).
2. **`estimator_influence_audit.R`** — the fairness audit, self-contained in R.
   For **both** acoustic spaces (2-D PB52 and 13-D SBCSAE MFCC) it correlates
   each separation metric against a panel of overlap yardsticks that do *not*
   share its machinery, and writes `estimator_influence_correspondence.csv`
   (dataset × metric × yardstick, |Spearman| and |Pearson| with bootstrap CIs).
   The neutral yardsticks (kNN density, MVN Monte-Carlo, PCA-2D grid / convex
   hull) are implemented in **`overlap_estimators.R`**, an R port of the poster's
   `reproduce_poster_figures.py`, validated to reproduce its values (rho > 0.98).
3. **`make_figures.R`** — renders PNG (300 dpi) + PDF to `outputs/figures/`:
   - `fig1_dimensionality_slopegraph` — correspondence with KDE overlap, 2-D vs 13-D
   - `fig2_mfcc_correspondence_panels` — per-metric scatter vs 1 − KDE overlap (13-D)
   - `fig3_estimator_robustness` — Monte-Carlo vs legacy JSD (the 2.1.0 fix preserves rank)
   - `fig4_estimator_influence_audit` — **centerpiece**: which metric matches
     overlap best depends on the overlap yardstick, in both dimensions

Each block/figure is skipped with a message if its input is missing.

## Inputs (place under `manuscript/data/`)

| file | source |
| --- | --- |
| `pb52.csv` | Peterson–Barney 1952 vowels (columns `vowel`, `f1`, `f2`, …); equivalently `phonTools::pb52` |
| `rerun_mfcc_sampled_tokens_20260624.csv` | committed on `labphon_2026` under `analysis/labphon_2026/data/` |

Everything else is computed from these tokens with phontrast 2.1.0 and
`overlap_estimators.R`; the audit needs no committed overlap values.

## Framing notes for the manuscript

- **The thesis is the package, not a single winning metric.** `phontrast`
  computes JSD/distance, Pillai, Bhattacharyya distance/affinity, Mahalanobis,
  and percent overlap side by side; the manuscript's contribution is enabling
  rigorous *comparative* analysis of acoustic contrast.
- **The estimator-influence audit is the utility demonstration (fig4), and it is
  a fairness check, not a bake-off.** JSD is estimated with KDE and so is
  phontrast's percent-overlap, so "JSD best matches KDE overlap" is partly
  circular. Correlating each metric against neutral overlap yardsticks shows the
  structure: **each metric matches best the yardstick built from its own
  assumptions** (JSD ↔ KDE, Pillai/Bhattacharyya ↔ MVN Monte-Carlo). On
  assumption-neutral references (kNN density, PCA) the metrics are comparable
  (mean |Spearman| ≈ 0.91 for all three) — no single metric is universally best.
- **The honest case for JSD is that it makes the fewest distributional
  assumptions of the separation metrics** (verified against the implementations):
  - JSD / Jensen–Shannon distance (`jsd_kde_nd`) are **nonparametric** — they
    read the whole distribution (shape, variance, multimodality) off a KDE and
    assume no parametric family.
  - Pillai (`pillai_overlap`, a MANOVA), Bhattacharyya (`bhattacharyya_mvnorm`,
    the Gaussian closed form), and Mahalanobis (`.mahalanobis_distance`, mean
    separation scaled by pooled covariance) all rest on **multivariate normality
    and/or a common-covariance assumption**, and use only first/second moments.
  - This connects directly to the audit: Pillai and Bhattacharyya look best
    against the **MVN Monte-Carlo** yardstick *because* they share its
    multivariate-normal model — an assumption real, often skewed or bimodal
    vowel clouds routinely violate. JSD earns its correspondence without that
    assumption, so it neither gets an MVN reference's free boost nor its penalty
    when the data are non-Gaussian.
  - Honest caveat: JSD is not assumption-*free*. It trades parametric
    assumptions for density estimation — a bandwidth choice and the curse of
    dimensionality — which is why it is the KDE-coupled metric in fig4.

  Net: recommend JSD on **robustness / minimal assumptions** (plus bounded,
  symmetric, dimension-general), not on a correspondence-ranking win.
- **Estimator choice:** `mc` is primary (consistent estimate of the continuous
  JSD; 2.1.0 default), `legacy` reported as a robustness column. The manuscript's
  vowel categories are large (PB52 E/I: 152 + 152 tokens; each SBCSAE vowel ≈
  500), where `mc` does not floor and its finite-sample null bias is negligible.
