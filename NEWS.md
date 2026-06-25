# phonJSD 0.5.1

## LabPhon 2026

- Added a reproducible LabPhon 2026 poster analysis bundle under
  `analysis/labphon_2026/`.
- Added the final poster PDF, processed analysis inputs, regenerated figure
  scripts, and audited no-MMO overlap-estimator outputs.
- Updated DESCRIPTION metadata so built-package checks materialize package
  author and maintainer fields correctly.

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
