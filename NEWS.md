# phontrast 2.5.1

## `rank_contrasts()` now reproduces the published calibration

**`rank_contrasts()` results change relative to 2.5.0.** The ranking protocol
it implements (Berry, "Estimand or estimator? Comparing vowel overlap measures
against a known ground truth", JASA, under review, Sec. VII A) was calibrated on
phontrast's `method = "legacy"` kernel estimator. Its sample-size floors
(`protocol_floors()`), its 0.25 margin, its 0.99 ceiling, and its bandwidth
check all hold for that estimator: the self-normalized sample-point
construction that evaluates both categories' densities at the pooled tokens,
with no leave-one-out correction. 2.5.0 computed `sqrt_jsd`, `shared_mass`, and
the bandwidth check with the Monte-Carlo estimator (`method = "mc"`, partial
leave-one-out) and a per-category Scott bandwidth instead, so its numbers were
not on the scale the thresholds were set on. On the paper's 45
Peterson--Barney F1 x F2 vowel pairs, 2.5.0 put 20 pairs at the ceiling instead
of 23, set no pair aside instead of ɝ--u (rank shift 0.080 against the paper's
0.273), and read `sqrt_jsd` up to 0.11 higher on overlapping pairs (ɑ--ɔ: 0.730
against the published 0.627). 2.5.1 reproduces the paper's per-pair values to
floating-point precision: 45 pairs, 23 at the ceiling, 22 ranked, none flagged,
ɝ--u set aside. Pillai is unchanged.

- `recommended_estimator()` returns three new elements: `method = "legacy"`,
  `eval_on = "pooled"`, and `bracket_bw = "scott.pooled"`, at every
  dimensionality. `loo` is now `FALSE` throughout, because the study applied no
  leave-one-out correction (`legacy` has none); the earlier `TRUE` described a
  correction that the calibrated estimator never applied.
- `rank_contrasts()` routes the reported `sqrt_jsd` and `shared_mass`, and the
  halved and doubled `scott.diag` bracket of the bandwidth check, through
  `estimator$method`. Under `bracket_bw = "scott.pooled"` the bracket uses one
  diagonal Scott bandwidth selected on the pooled tokens of the pair for both
  categories, as the study did; `bracket_bw = "scott.diag"` selects it per
  category, as 2.5.0 did. The `estimator` checker validates the new elements;
  a list without them, such as one built by hand for 2.5.0, takes the
  calibrated values. To recompute a 2.5.0 ranking, pass
  `modifyList(recommended_estimator(d), list(method = "mc", loo = TRUE, bracket_bw = "scott.diag"))`.
- `print()` on a ranking names the estimator method and says when the bracket
  is pooled; `inspect_contrast()` draws and annotates its panels with the
  ranking's own estimator and bracket bandwidth, so the x0.5 and x2 panels
  repeat `sqrt_jsd_half` and `sqrt_jsd_double`.
- On the bundled `vowel_cohort`, the bandwidth check now sets aside spk09 (the
  planted disagreement, still flagged) and spk06; the vignette says why.
- A regression test reproduces the paper's Peterson--Barney outcome and
  per-pair values from a vendored fixture (`tests/testthat/fixtures/`).

## Documentation

- The help for `method = "legacy"` in `jsd_kde_nd()` and its siblings
  (`estimate_jsd()`, `percent_overlap_kde()`, `estimate_overlap()`,
  `phontrast()`, `compare_overlap_metrics()`, `jsd_summary()`,
  `global_boot_jsd()`) no longer calls it a pre-1.2.0 index to use only to
  reproduce phonJSD 1.0.0. It now says that `legacy` is the estimator the
  published protocol was validated on, and that no leave-one-out correction
  applies under it.

## Other changes

- The metric annotation boxes drawn by `plot_contrast()` and
  `inspect_contrast()` no longer use `geom_label(label.size =)`, which recent
  ggplot2 releases deprecate; they pass the `linewidth` aesthetic when the
  installed ggplot2 has it and `label.size` otherwise. Plots are unchanged.
- The citation test no longer depends on where the text rendering wraps its
  lines, which a four-component development version number could trip.

# phontrast 2.5.0

## The ranking protocol in one call

This release implements the measurement protocol recommended by the simulation
study "Estimand or estimator? Comparing vowel overlap measures against a known
ground truth" (Berry, under review, Sec. VII.A), so that its three steps and
their conditions can be run and reported without hand-assembling them.

- `rank_contrasts(data, features, category_col, group_col)` computes, per
  speaker, Jensen-Shannon distance (`sqrt_jsd`) and shared probability mass
  from one shared kernel density estimate and the Pillai trace on the same
  tokens (step 1); ranks the speakers on the percentile-rank scale
  `(average rank - 1/2) / k`, ties averaged, after setting measurements at the
  `sqrt(JSD)` ceiling (default 0.99) apart (step 2); and flags speakers whose
  Pillai and `sqrt(JSD)` percentile ranks differ by `margin` (default 0.25 of
  the ordering) or more (step 3).
- The study's conditions are applied and reported rather than left to the
  user: sample-size floors per speaker from the smaller category's count
  (`rank_licensed`, `flag_licensed`, `rank_basis`; ranking licensed from 50
  tokens per category at d = 2, 200 at d = 3-4, 500 at d = 5-8; the flag
  readable from 100 at d = 2 and 200 at d = 3-4, undefined when `d >= 2n`,
  withheld from d = 5), a warning below eight speakers, and the bandwidth
  check (`sqrt(JSD)` recomputed at half and twice the diagonal Scott
  bandwidth; `bw_shift`, `sign_change`, `set_aside`).
- The result is a `phontrast_ranking` tibble with a `print()` method that
  summarizes the protocol, the estimator, the floors, the flagged and
  set-aside speakers, and the ceiling. The cleaned tokens travel with it in
  `attr(x, "protocol")`.
- `plot_rank_agreement()` (also `plot()` / `ggplot2::autoplot()` on the
  ranking) draws Pillai percentile rank against `sqrt(JSD)` percentile rank
  with the identity line and the `+/- margin` inspection band; flagged speakers
  are coloured and labelled, set-aside speakers crossed, unreadable flags
  hollow, and ceiling speakers listed in the caption.
- `inspect_contrast(ranking, group)` redraws one speaker with
  `plot_contrast()`'s distribution-aware layers at half, the selected, and
  twice the diagonal Scott bandwidth, annotating each panel with `sqrt(JSD)`
  and shared mass at that bandwidth plus the speaker's Pillai, and restating
  the reported ranks and the bandwidth-check outcome in the subtitle.
- Helpers: `percentile_rank()` (the protocol's rank scale, with an `exclude`
  mask), `protocol_floors(d)` (the licensing floors by dimensionality), and
  `recommended_estimator(d)` (the kernel estimator settings the study used at
  each dimensionality: plug-in / `ks` at d <= 4, diagonal Scott / `fast_diag`
  with 200 evaluation tokens from d = 5, no leave-one-out from d = 14).
- A bundled, simulated twelve-speaker cohort, `vowel_cohort`, with a planted
  Pillai / `sqrt(JSD)` disagreement, a ceiling speaker, and two under-sampled
  speakers; and a vignette, "Ranking speakers by Jensen-Shannon distance and
  checking Pillai agreement", that walks through the protocol on it.

## Bandwidth multiplier

- New `bw_scale` argument on the kernel path (`jsd_kde_nd()`,
  `percent_overlap_kde()`, `estimate_jsd()`, `estimate_overlap()`,
  `phontrast()`, `plot_contrast()`): a positive multiplier on the selected
  bandwidth on the standard-deviation scale (`h * bw_scale`; bandwidth
  matrices `H * bw_scale^2`). `0.5` and `2` give the halved and doubled
  bandwidths of the smoothing-sensitivity check. `plot_contrast()` draws its
  regions, overlap shading, and annotations at the scaled bandwidth and
  records it in the caption. Ignored under `density = "mvnorm"`.

## The kernel family on one shared density

- `phontrast()` now reads Jensen-Shannon divergence and distance, proportional
  overlap, total variation, and the matched-kernel Bhattacharyya affinity /
  distance and Hellinger distance off a single density estimate per
  comparison, instead of separate kernel passes for JSD and overlap. Values of
  the existing columns are unchanged; bootstraps do half the kernel work.
- New opt-in metrics for `phontrast(metrics = ...)`: `"tv"` (total variation,
  `total_variation = 1 - percent_overlap`); `"bhattacharyya_kde"`
  (`bhatt_kde_dist`, `bhatt_kde_affinity`, `hellinger`), the same quantity as
  the closed-form Gaussian `"bhattacharyya"` columns under a matched kernel
  estimator, with the same partial leave-one-out self-density as JSD, so
  measure and estimator can be told apart; and `"euclidean"`
  (`euclidean_dist`), the distance between the category means after dividing
  each feature by the pooled two-category standard deviation (bounded by
  `2 * sqrt(d)` for equal category sizes). The default metric set and the
  column order of earlier releases are unchanged; long output and bootstraps
  cover the new columns.

## Citation

- `citation("phontrast")` now lists the study behind the protocol as a second
  entry (manuscript under review, with the OSF replication package), after the
  software citation.

# phontrast 2.4.1

## Cross-BLAS robustness of the proportion-standardized Pillai guard

- Fixed the test failure reported by CRAN's tests-MKL additional check. The
  nonsingularity guard for proportion-standardized Pillai estimates decided
  rank deficiency by whether `chol()` of the within-class error SSCP threw an
  error, which is BLAS-dependent on (near-)singular input: reference LAPACK
  errors where Intel MKL can return a tiny positive pivot, letting degenerate
  designs slip past the guard and fail later inside `summary.manova()` with a
  less specific message. The guard now decides rank by R's tolerance-based QR
  of the residual matrix -- the same criterion `summary.manova()` applies --
  so degenerate designs raise the documented "nonsingular" error on every
  BLAS build. Estimates on well-posed designs are unchanged.
- Added a near-collinear regression test that pins this behavior.

# phontrast 2.4.0

## Version metadata

- Corrected the package version in `DESCRIPTION` to 2.4.0. The original
  release tags contained the 2.4.0 changes but retained the 2.3.1 package
  metadata.
- Replaced the citation that was frozen at v2.0.0 with version-aware package
  citation metadata and the version-specific Zenodo DOI. Added `CITATION.cff`
  so GitHub and future Zenodo deposits use the package's declared authorship.

## Proportion-standardized Pillai estimates

- Added the opt-in `proportion_standardized = TRUE` mode to
  `pillai_overlap()`. It appends the plug-in and unbiased squared Mahalanobis
  separation estimates and Becker's proportion-standardized Pillai score while
  leaving the default two-field return and raw `pillai` value unchanged.
- Negative unbiased separation estimates follow an explicitly labelled
  fallback path: `pillai_eq` is `NA`, `pillai_eq_fallback` is `TRUE`, and the
  multiplicatively corrected first term appears only as `d2_fallback` beside
  its closed-form upward-bias term.
- Added definedness and fragility diagnostics for the two-category,
  no-covariate estimator chain, along with regression vectors, base-R MANOVA
  comparisons, closed-form anchors, and null Monte-Carlo acceptance tests.

# phontrast 2.3.1

## CRAN release

- Published on CRAN for the first time on 2026-08-09. The stable release is
  available with `install.packages("phontrast")`.

## Proportion-standardized Pillai estimates

- Added the opt-in `proportion_standardized = TRUE` mode to
  `pillai_overlap()`. It appends the plug-in and unbiased squared Mahalanobis
  separation estimates and Becker's proportion-standardized Pillai score while
  leaving the default two-field return and raw `pillai` value unchanged.
- Negative unbiased separation estimates follow an explicitly labelled
  fallback path: `pillai_eq` is `NA`, `pillai_eq_fallback` is `TRUE`, and the
  multiplicatively corrected first term appears only as `d2_fallback` beside
  its closed-form upward-bias term.
- Added definedness and fragility diagnostics for the two-category,
  no-covariate estimator chain, along with regression vectors, base-R MANOVA
  comparisons, closed-form anchors, and null Monte-Carlo acceptance tests.

## CRAN resubmission fixes

- Replaced the non-running `extract_mfcc()` example with a fast, executable
  example that creates and removes a short WAV file in the R session's
  temporary directory. The optional `tuneR` dependency is guarded with
  `requireNamespace()` as recommended for packages in `Suggests`.
- Unwrapped the short `phontrast()` and `hier_boot_jsd_model()` bootstrap
  examples after confirming that they run well under five seconds. Reduced the
  illustrative `estimate_jsd()` bootstrap count so that example is also under
  five seconds and can be unwrapped. No examples now use `\dontrun{}` or
  `\donttest{}`.
- Removed all direct access to `.GlobalEnv`. Seeded KDE subsampling and
  multivariate-normal Monte Carlo now use a private deterministic generator,
  retaining reproducibility without reading, writing, or replacing the user's
  `.Random.seed`. Unseeded stochastic calls retain their previous behavior.

# phontrast 2.3.0

## Distribution-aware, accountable plotting

- **New flagship visualization `plot_contrast()`.** Draws the *same density
  model the distributional metrics are computed from*: under
  `density = "kde"` it shows highest-density regions of each category's kernel
  density estimate (same bandwidth selection and ks evaluation as the
  metrics); under `density = "mvnorm"` it shows coverage ellipses of the
  fitted Gaussians used by the parametric backend. The pointwise minimum of
  the two densities -- the mass the proportional-overlap metric integrates --
  is shaded directly (a ribbon in 1-D, a soft raster in 2-D, normalized across
  panels so fainter panels genuinely overlap less).
- **Accountability by default.** Panels are annotated with the Jensen-Shannon
  divergence and proportional overlap computed by `phontrast()` under the
  plotted density model (optionally with bootstrap intervals via `n_boot`);
  the caption records the estimator configuration (backend, bandwidth,
  `mc_n`, seed, n); and the full annotation table is attached to the plot as
  `attr(p, "contrast_metrics")`.
- **Results plot themselves.** `phontrast()` output now carries class
  `"phontrast_contrast"`, so `plot()` and `ggplot2::autoplot()` draw the
  metric comparison directly via `plot_overlap_metrics()`.
- **A shared visual identity.** New exported `theme_phontrast()`,
  `scale_colour_phontrast()` / `scale_color_phontrast()` /
  `scale_fill_phontrast()`, and `phontrast_palette()` (a colorblind-safe
  Okabe-Ito palette ordered so the leading pair maximizes contrast for
  two-category plots). All plotting functions -- including the existing
  `plot_overlap_metrics()`, `plot_category_space()`, and
  `plot_category_pca()`, whose arguments and behavior are otherwise
  unchanged -- now share this theme and palette.
- Plotting remains optional: `ggplot2` stays in Suggests.

# phontrast 2.2.0

## Pluggable density backend for the distributional metrics

- **Added a `density` argument to decouple the density *estimator* from the
  *metric*.** The distributional metrics (Jensen-Shannon divergence and
  proportional overlap) previously always used kernel density estimation. They
  now accept `density = "kde"` (the default, unchanged) or `density = "mvnorm"`,
  which fits one multivariate normal per category and estimates the metric
  between the two Gaussians. This lets the estimator behind JSD and overlap be
  matched to the same multivariate-normal assumptions the Pillai, Bhattacharyya,
  and Mahalanobis metrics already make, and makes the metric x estimator
  interaction a controlled choice rather than hard-wired to KDE.
- Jensen-Shannon divergence between two Gaussians has no closed form (the
  mixture is a Gaussian mixture), so the `"mvnorm"` backend estimates it by
  **fresh-sample Monte-Carlo**: it draws `mc_n` points (default `10000`) from
  each fitted Gaussian and averages the log density ratio. The estimand is the
  JSD / overlap between the fitted Gaussians; `eval_seed` makes the draw
  reproducible without disturbing the caller's random-number stream. The
  Gaussian fit has no self-kernel, so the KDE-specific leave-one-out correction
  does not apply.
- `density` (and `mc_n`) are threaded through `phontrast()`,
  `compare_overlap_metrics()`, `estimate_jsd()`, `estimate_overlap()`,
  `jsd_summary()`, `global_boot_jsd()`, `jsd_kde_nd()`, and
  `percent_overlap_kde()`. In `phontrast()` the argument affects only the
  Jensen-Shannon and overlap columns; the Pillai, Bhattacharyya, and
  Mahalanobis columns are parametric by construction and are unchanged.
- The multivariate-normal log-density and sampler are implemented in base R via
  a Cholesky solve, adding no new package dependency. Rank-deficient category
  covariances are regularized with a small ridge before factorization.
- Fully backward compatible: the default remains `density = "kde"` and all
  existing results are unchanged.

# phontrast 2.1.0

## Monte-Carlo JSD no longer floors small divergences to exactly 0

- **Fixed `method = "mc"` (the default KDE JSD estimator) flooring small but
  real divergences to exactly 0.** The full leave-one-out correction could
  collapse a category's self-density at isolated points, driving the raw
  plug-in mean negative, which the final clamp then floored to 0 -- while
  `method = "legacy"` still reported a nonzero contrast. The estimator now uses
  a *partial* leave-one-out correction that removes a sample-size-scaled
  fraction `n / (n + 20)` of each point's own kernel: half at 20 tokens per
  category (the `min_tokens` default), approaching the full correction as the
  category grows. The corrected density stays strictly positive, so
  near-merged categories yield small positive estimates instead of exact 0.
- Validated against grid-integrated JSD of the same KDEs: on 44 real
  speaker contrasts (20 tokens per phase, 2-D formant space) the previous
  estimator returned exactly 0 for 24 speakers; the corrected one returns 0
  for none, tracks the grid reference more closely than either the previous
  default or `method = "legacy"`, and still passes the package's
  grid-calibration test at n = 200 within the original tolerance.
- **This changes `method = "mc"` results relative to 2.0.x**, most visibly for
  small divergences and small samples (estimates that were floored at 0 become
  small positive values; others typically shift upward slightly).
  `method = "legacy"` is unchanged, and `loo = FALSE` is unchanged.

## CRAN preparation

- Replaced the relative `ROADMAP.md` link in the README with a plain reference,
  since `ROADMAP.md` is excluded from the built package; this resolves the
  "invalid file URI" flagged by the CRAN incoming checks.

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
