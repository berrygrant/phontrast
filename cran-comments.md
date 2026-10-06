## Submission

phontrast 2.5.0 is a feature release. It adds a one-call implementation of a
measurement protocol for ranking speakers' vowel contrasts by Jensen-Shannon
distance and checking the ranking against the Pillai trace
(`rank_contrasts()`, with sample-size licensing, a bandwidth-sensitivity
check, a rank-agreement plot, and an inspection plot), a bandwidth multiplier
across the kernel-density path, three opt-in metrics for `phontrast()`, a
small simulated example dataset (`vowel_cohort`, 2200 rows), and a vignette
that walks through the protocol. The kernel-family metrics that `phontrast()`
already reported are now computed from one shared density estimate per
comparison; their values are unchanged. There are no API removals.

All new examples run in under 2.5 seconds each on the local machine (the
heaviest, `inspect_contrast()`, in 2.2 s); the release keeps the example
budget that was corrected for 2.4.0.

## Test environments

- Local: Ubuntu 24.04, R 4.3.3 with reference BLAS/LAPACK (offline build
  environment; the network-dependent CRAN incoming checks could not run there)
- R-hub v2 (GitHub Actions): Ubuntu R-devel, Windows R-devel, macOS R-release,
  and the `mkl` (Intel MKL) container
- GitHub Actions: Ubuntu latest, R release, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

The sole local note is environmental: the suggested package `tuneR` is not
installable in the offline build container, so it was unavailable for
checking there. It is available on CRAN's machines.

The local PDF-manual check was skipped (no TeX installation); all Rd checks,
examples, tests (649 expectations), and the three vignettes completed
successfully.
