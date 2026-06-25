# phonJSD 0.5.1

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
