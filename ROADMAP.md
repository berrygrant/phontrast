# phontrast roadmap

phontrast grew out of `phonJSD`, a package focused on Jensen-Shannon
divergence, and was reoriented into a general toolkit for computing and
comparing multiple phonological category **contrast and separation metrics**.
This document records where it is and where it is headed. Last revised for
2.5.0.

## Shipped — the 2.x series

The 2.0.0 plan ("P1 — full architectural redesign") set out eight steps. Most
of what users were waiting for has shipped as additive releases; the
architectural steps that remain are listed under *3.0.0* below.

| Release | What shipped | Plan item |
| --- | --- | --- |
| 2.0.0 | Renamed `phonJSD` → `phontrast`; `phontrast()` as the unified entry point for any subset of the metrics, globally or by group, with bootstrap intervals; `compare_overlap_metrics()` deprecated; corrected KDE estimator. | — |
| 2.1.0 | Monte-Carlo JSD with the sample-size-scaled partial leave-one-out correction, so small real divergences are no longer floored to exactly 0. | — |
| 2.2.0 | Pluggable density backend: `density = c("kde", "mvnorm")` on `phontrast()`, `estimate_jsd()`, `estimate_overlap()`, `jsd_summary()`, `global_boot_jsd()`, `jsd_kde_nd()`, `percent_overlap_kde()`; the Gaussian backend estimates JSD and overlap between the fitted Gaussians by fresh-sample Monte-Carlo (`mc_n`, `eval_seed`). | 3 (KDE and MVN; GMM still open) |
| 2.3.0 | Distribution-aware, accountable plotting: `plot_contrast()` draws the same density model the metrics use and annotates panels with their values; `plot()` / `autoplot()` on `phontrast()` results; the Okabe-Ito `theme_phontrast()` family; `plot_overlap_metrics()`, `plot_category_space()`, `plot_category_pca()` retrofitted. | *Later: richer visualisation* |
| 2.3.1 | First CRAN release (2026-08-09). Opt-in proportion-standardized Pillai (`pillai_overlap(proportion_standardized = TRUE)`). | 8 |
| 2.4.0 / 2.4.1 | Version-aware `inst/CITATION`, `CITATION.cff`, Zenodo DOIs per release; cross-BLAS robustness of the proportion-standardized Pillai guard (CRAN tests-MKL). | — |
| 2.5.0 | The measurement protocol of the JASA simulation study (Sec. VII.A) in one call: `rank_contrasts()` with `percentile_rank()`, `protocol_floors()`, `recommended_estimator()`, `plot_rank_agreement()` (`plot()` on the ranking), and `inspect_contrast()`; `bw_scale` across the kernel path; the kernel family (JSD, overlap, total variation, kernel Bhattacharyya / Hellinger) scored on one shared density per comparison; opt-in `"tv"`, `"bhattacharyya_kde"`, `"euclidean"` metrics; bundled `vowel_cohort`; protocol vignette. | *Later: additional metrics* (in part) |

Two further plan items are in place in substance, if not in the form the
plan described: the bootstrap resamples every requested metric with uniform
`*_mean`, `*_sd`, `*_ci_lower`, `*_ci_upper` columns (item 4), and long
output carries `orientation`, `separation_value`, and `separation_rank`
(item 5). The metric list in `phontrast()` is an implicit registry (items 1
and 2, partially): one internal table names each metric's columns and
orientation, but adding a metric still touches several places.

## Next — 2.6.x (additive)

- **The rest of the study's measure set.** SOAM (Wassink 2006; 2-SD ellipses
  on the covariance principal axes, defined at two and three dimensions) and
  APP (Morrison 2008; per-category QDA trained on fresh draws from the fitted
  Gaussians) as opt-in `phontrast()` metrics, so the paper's full seven-measure
  table can be reproduced in-package. APP also covers the "classifier-based
  separability" idea from the earlier roadmap.
- **Protocol follow-ups.** A report helper that emits the "what to report"
  block of the protocol vignette as text or a table; several contrasts per
  speaker in one `rank_contrasts()` call; `summary()` for `phontrast_ranking`.
- **Overlap estimator consistency.** `percent_overlap` evaluates the raw
  self-density where JSD and the kernel Bhattacharyya use the partial
  leave-one-out correction. Aligning the three is a results-changing step
  for `percent_overlap` and will be flagged in NEWS if taken.
- **Gaussian-mixture density backend** (`density = "gmm"`): a larger step
  that adds an EM dependency and component selection.
- Energy distance as a further distribution-free metric.

## 3.0.0 — the architectural steps

These change internals and remove shims, so they wait for a major version.

1. **Formal metric registry and per-metric contract** (plan items 1, 2, 7).
   Register each metric with id, label, orientation, theoretical range,
   bootstrap support, and modelling assumptions; `phontrast()` dispatches
   through the registry; adding a metric means registering one
   `estimate(data, features, group, ...) -> scalar` function; document how
   users register their own.
2. **Consolidate wrappers** (plan item 6). Refactor `speaker_*`,
   `estimate_*`, `jsd_summary()`, and `hier_boot_jsd_model()` into thin shims
   over the unified core, deprecate gradually, and remove in 3.0.0.

## Maintenance

- Per release: the Zenodo DOI into `inst/CITATION` and `test-citation.R`,
  the README install pin, and the `CITATION.cff` release date.
- Replace the "under review" entry in `inst/CITATION` with the published
  reference once the JASA paper is accepted.

Contributions and suggestions are welcome via
<https://github.com/berrygrant/phontrast/issues>.
