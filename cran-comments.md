## Submission

This is a new submission. `phontrast` is a rename and reframing of the
previously GitHub-only package `phonJSD` (released through v1.2.0); it has not
been on CRAN before. Development, testing, and documentation were assisted by AI
tools, as disclosed in the README; all design, analysis, and release decisions
were made by the maintainer.

## Test environments

<!-- Fill in with your actual runs before submitting. -->

- Local: <your OS>, R <version>
- win-builder: R-devel and R-release (`devtools::check_win_devel()`,
  `check_win_release()`)
- macOS builder: <mac.r-project.org result>
- R-hub: <platforms>

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
