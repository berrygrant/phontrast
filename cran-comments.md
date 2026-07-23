## Resubmission

This is a resubmission. Compared to the previously submitted version (2.0.2):

- I removed a relative-path link to `ROADMAP.md` from the README (that file is
  intentionally not shipped in the package), which resolves the "invalid file
  URI" reported by the CRAN incoming checks.
- The version was bumped to 2.1.0 because the default Monte-Carlo
  Jensen-Shannon divergence estimator was corrected: its leave-one-out bias
  correction could floor small but real divergences to exactly 0, and now uses
  a sample-size-scaled partial correction instead (documented in NEWS.md and
  covered by a regression test).

## Submission

This is a new submission (first time on CRAN). `phontrast` is a rename and
reframing of the previously GitHub-only package `phonJSD` (released through
v1.2.0). Development, testing, and documentation were assisted by AI tools, as
disclosed in the README; all design, analysis, and release decisions were made
by the maintainer.

## Test environments

- Local: macOS 26.5 (aarch64-apple-darwin), R 4.5.3
- win-builder: R-devel (`devtools::check_win_devel()`)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Notes for the reviewer

* The spell-check flags domain terminology in the Description and
  documentation (e.g. "Bhattacharyya", "Mahalanobis", "Pillai", "MFCCs",
  "sociophonetics", "phonological"). These are standard terms from phonetics
  and multivariate statistics; they are listed in `inst/WORDLIST`.
* The `Description` cites Lin (1991) <doi:10.1109/18.61115> for the
  Jensen-Shannon divergence.
* `ggplot2`, `mgcv`, and `tuneR` are used only conditionally
  (`requireNamespace()`), so they are in Suggests rather than Imports.
