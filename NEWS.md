# phontrast 2.0.2

## CRAN preparation

- Kept non-package content out of the built tarball: the (untracked) local
  `OSF/` reproducibility-data directory and the top-level `ROADMAP.md` are now
  listed in `.Rbuildignore`, resolving the non-portable-paths and
  non-standard-top-level-files NOTEs from `R CMD check`. Added "mel" to
  `inst/WORDLIST`.
- Trimmed the `estimate_jsd()` example to a fast point estimate and moved the
  bootstrap demonstrations into `\donttest{}`, keeping every example under
  CRAN's execution-time limit.

# phontrast 2.0.1

## CRAN preparation

- Prepared the package for CRAN submission: expanded acronyms and added a
  Jensen-Shannon divergence reference (Lin, 1991) to the `Description`,
  normalized non-ASCII characters in the R sources, added `inst/WORDLIST` and
  `cran-comments.md`, and removed the AI assistant from `Authors@R` (the AI-use
  acknowledgment remains in the README). No user-facing code changes.

# phontrast 2.0.0

## Package renamed: phonJSD is now phontrast

- **The package has been renamed from `phonJSD` to `phontrast`** and reoriented
  around comparing *multiple* category contrast and separation metrics rather
  than Jensen-Shannon divergence alone. Update your code from `library(phonJSD)`
  to `library(phontrast)`. Function names are unchanged except as noted below,
  and no metric estimates change relative to 1.2.0 -- this release is a rename
  and API reframe, not a numerical change.
- The GitHub repository moved to <https://github.com/berrygrant/phontrast> (old
  links redirect), and the Zenodo concept DOI (10.5281/zenodo.20816585) is
  unchanged, so existing citations continue to resolve.

## New unified entry point: `phontrast()`

- Added **`phontrast()`**, the package's headline function: compute and compare
  any subset of the contrast metrics -- Jensen-Shannon divergence and distance,
  Pillai-Bartlett trace, Bhattacharyya distance and affinity, Mahalanobis
  distance, and proportional overlap -- for a two-category contrast in one call,
  globally or by group, wide or tidy long, with optional bootstrap intervals.
- The new `metrics` argument selects which metrics to compute (default: all),
  e.g. `phontrast(data, features, "vowel", metrics = c("jsd", "pillai"))`.
- **`compare_overlap_metrics()` is deprecated** in favor of `phontrast()`. It
  still works (it calls `phontrast()` with `output = "wide"`) and emits a
  deprecation message; it will be removed in a future release.

## Roadmap

- The next priority (P1) is a full architectural redesign around a metric
  registry with a uniform per-metric interface, a generalized bootstrap that
  works for any metric, and first-class metric orientation. See `ROADMAP.md`.

# phonJSD 1.2.0

## Corrected KDE estimator (changes results)

- **KDE-based JSD and percent overlap now use a consistent Monte-Carlo plug-in
  estimator by default (`method = "mc"`).** Each category's KDE is evaluated at
  that category's own observations and the true log density ratio against the
  mixture is averaged (with a leave-one-out bias correction, `loo = TRUE`). This
  estimates the continuous Jensen-Shannon divergence in any dimension, replacing
  the previous self-normalized sample-point index, which was a bounded relative
  separation measure rather than the JSD integral and depended on `eval_on`.
- **This changes the numbers relative to phonJSD 1.0.0.** To reproduce 1.0.0
  results exactly, pass `method = "legacy"` to `jsd_kde_nd()`,
  `percent_overlap_kde()`, `estimate_jsd()`, `estimate_overlap()`,
  `jsd_summary()`, `global_boot_jsd()`, or `compare_overlap_metrics()`. The
  `eval_on` control applies to `method = "legacy"` only.
- The `fast_diag` engine now evaluates true (normalized) densities, matching
  `ks::kde()` to machine precision for diagonal bandwidths.

## API

- Aligned defaults across the estimation API: `min_tokens = 20` and
  `n_boot = 1000` everywhere (previously `estimate_jsd()` defaulted to
  `min_tokens = 5`, and `jsd_summary()`/`boot_jsd()`/`hier_boot_jsd_model()` to
  `min_tokens = 30` / `n_boot = 300`).
- Standardized the `estimate_*`/`global_*` wrappers to return tibbles uniformly.
- Documented the `pillai_p_value` (wide) and `p_value` (long) columns returned
  by `compare_overlap_metrics()`.

## Metrics

- Added opt-in high-dimensional KDE speed controls: `bw = "scott.diag"`,
  evaluation-point subsampling via `eval_n`/`eval_seed`, and
  `engine = "fast_diag"` for chunked diagonal-Gaussian KDE evaluation.
- Added `engine = "fast_diagonal"` as an alias for `engine = "fast_diag"`.
- Extended the same KDE controls across JSD, percent overlap, and
  `compare_overlap_metrics()` so separation and overlap estimates use aligned
  density-estimation settings.
- Exposed the same KDE controls in lower-level JSD wrappers including
  `speaker_jsd()`, `boot_jsd()`, `jsd_summary()`, and `global_boot_jsd()`.
- Allowed grouped metric wrappers to use multiple grouping columns via
  `group_col = c("Sex", "Style")`; grouped outputs retain a single labeled
  `group` column.
- Switched `extract_mfcc()` to `tuneR::melfcc()` and removed the stale
  `seewave::mfcc()` reference.

## Documentation

- Documented the high-dimensional fast KDE path in the README and generated
  function manuals.
- Updated the GitHub Actions checkout step to a current Node runtime action.

## Bug fixes

- Clamped `jsd()` output to the mathematical range `[0, 1]`, so floating-point
  rounding on near-identical categories can no longer yield `NaN`
  Jensen-Shannon distances (`est_distance = TRUE`) or abort
  `hier_boot_jsd_model()` through `prepare_jsd_beta()`.
- Fixed a latent `sample()` edge case in `hier_boot_jsd_model()` that could
  misdraw a single numeric group identifier.

## Robustness

- Grouped metric wrappers now keep a group as `NA` when its metric cannot be
  computed and emit a single summarizing warning, instead of the previous
  inconsistent behavior where `speaker_pillai()`/`speaker_bhatt()` silently
  dropped failed groups (so `estimate_pillai()` and `estimate_jsd()` could
  return different rows for the same data) while other metrics returned a silent
  `NA`.
- Added input validation at the main entry points: `data` must be a data frame,
  `category_col` a single column name, `features` a non-empty character vector,
  and `features` may not overlap with `category_col`/`group_col`.

## Package quality

- Removed the unused `LazyData` field (no `data/` directory), the unused `lme4`
  suggestion, and dead `dplyr` imports (`filter`, `n`, `ungroup`); added `URL`
  and `BugReports`; standardized the author name to "Grant M. Berry"; and scoped
  the CI workflow to the existing `main` branch.

# phonJSD 1.0.0

## Visualization

- Added ggplot2-backed `plot_overlap_metrics()` and `plot_category_space()`
  helpers for visualizing metric comparisons and one- or two-dimensional
  phonological category spaces.
- Added `plot_category_pca()` for two-dimensional PCA diagnostics of arbitrary
  multidimensional feature spaces.

## Metrics

- Added `compare_overlap_metrics()` to compute Pillai trace, Bhattacharyya
  distance and affinity, Jensen-Shannon divergence and distance, Mahalanobis
  distance, and percent overlap in one global or grouped comparison table.
- Added one-dimensional KDE support for JSD and percent-overlap estimates.
- Aligned KDE bandwidth/evaluation controls across JSD and percent overlap.
- Made grouped JSD bootstrap intervals respect `conf_level`.
- Improved empty grouped outputs and small-sample diagnostics across metric
  wrappers.
- Made two-category metric checks ignore unused factor levels after filtering,
  so filtered factor data such as PB52 `I/i` contrasts work without manually
  calling `droplevels()`.
- Clarified that percent-overlap outputs are 0--1 proportions, not 0--100
  percentages.

## Documentation

- Added quick-start and multidimensional-workflow vignettes.
- Added a README quick start that starts from a vowel token table and moves to
  metric output.
- Added metric-choice guidance and PB52 small-sample notes.
- Added preferred-path and metric-direction guidance to the package manual.
- Added runnable examples for the main metric comparison, JSD estimation,
  direct KDE JSD, beta-regression preparation, MFCC extraction, and
  hierarchical bootstrap modeling workflows.
- Standardized JSD bootstrap outputs with `conf_level`, `ci_lower`, and
  `ci_upper` columns while retaining `jsd_low` and `jsd_high` as aliases.
- Added optional bootstrapping to `compare_overlap_metrics()` with progress
  messages and metric-specific confidence intervals.

## Maintenance

- Aligned DESCRIPTION metadata with the v1.0.0 README/DOI release.
- Added explicit Author and Maintainer fields for source-package checks on
  current R tooling.
- Removed old draft analysis scripts from the main package branch; the final
  LabPhon 2026 poster and reproducibility bundle live on `labphon_2026`.

# phonJSD 0.5.0

## Package Quality

- Added testthat coverage for the discrete JSD core, grouped classical metrics,
  KDE-based JSD wrappers, and grouped bootstrap summaries.
- Added `.Rbuildignore` entries so local analysis, load-test, and manuscript
  artifacts are excluded from package builds.
- Fixed the MIT license stub expected by R package tooling.

## Metrics

- Centralized column validation, complete-case filtering, and two-category
  checks across metric functions.
- Corrected discrete KL/JSD handling of zero-probability events.
- Made grouped bootstrap JSD report the number of successful bootstrap
  replicates.
- Replaced bootstrap `replicate()` control-flow edge cases with explicit
  `vapply()` iteration.
- Removed deprecated tidyselect usage in grouped JSD summaries.

## Documentation

- Updated generated Rd documentation for KL/JSD and bootstrap JSD outputs.
- Updated README release metadata and feature summary for v0.5.0.
