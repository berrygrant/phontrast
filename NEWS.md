# phonJSD 0.5.1

## Metrics

- Added `compare_overlap_metrics()` to compute Pillai trace, Bhattacharyya
  distance and affinity, Jensen-Shannon divergence and distance, Mahalanobis
  distance, and percent overlap in one global or grouped comparison table.
- Added one-dimensional KDE support for JSD and percent-overlap estimates.
- Aligned KDE bandwidth/evaluation controls across JSD and percent overlap.
- Made grouped JSD bootstrap intervals respect `conf_level`.
- Improved empty grouped outputs and small-sample diagnostics across metric
  wrappers.

## Documentation

- Added a README quick start that starts from a vowel token table and moves to
  metric output.
- Added preferred-path and metric-direction guidance to the package manual.
- Added runnable examples for the main metric comparison, JSD estimation,
  direct KDE JSD, beta-regression preparation, MFCC extraction, and
  hierarchical bootstrap modeling workflows.

## Maintenance

- Aligned DESCRIPTION metadata with the 0.5.1 README/DOI release.
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
