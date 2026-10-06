# Builds data/vowel_cohort.rda: a simulated twelve-speaker cohort for the
# ranking-protocol vignette and examples (rank_contrasts()). Run from the
# package root. The cohort is deterministic (fixed seed).
#
# spk01-spk08  Gaussian categories, centroid gap growing from 25 to 200 Hz,
#              100 tokens per vowel: a graded ordering both measures recover.
# spk09        Same centroids for both vowels, but a bimodal "eh" (two variants
#              190 Hz apart in F1, 440 Hz apart in F2): a mean-based measure
#              sees no contrast while the distributions barely overlap. This is
#              the planted Pillai / sqrt(JSD) disagreement the flag should catch.
# spk10        Fully separated categories: sqrt(JSD) at the ceiling.
# spk11        60 tokens per vowel: ranking licensed at two dimensions, flag not
#              readable (floor 100).
# spk12        40 tokens per vowel: below the two-dimensional rank floor (50).

set.seed(20260514)

make_speaker <- function(id, n, gap, sd1 = 55, sd2 = 120, f1_0 = 500, f2_0 = 1900) {
  data.frame(
    speaker = id,
    vowel = rep(c("ih", "eh"), each = n),
    f1 = c(rnorm(n, f1_0, sd1), rnorm(n, f1_0 + gap, sd1)),
    f2 = c(rnorm(n, f2_0, sd2), rnorm(n, f2_0 - gap, sd2)),
    stringsAsFactors = FALSE
  )
}

n_full <- 100
gaps <- c(25, 50, 75, 100, 125, 150, 175, 200)
graded <- Map(make_speaker, sprintf("spk%02d", 1:8), n_full, gaps)

side <- sample(c(-1, 1), n_full, replace = TRUE)
spk09 <- data.frame(
  speaker = "spk09",
  vowel = rep(c("ih", "eh"), each = n_full),
  f1 = c(rnorm(n_full, 500, 30), 500 + 95 * side + rnorm(n_full, 0, 30)),
  f2 = c(rnorm(n_full, 1900, 70), 1900 - 220 * side + rnorm(n_full, 0, 70)),
  stringsAsFactors = FALSE
)
spk10 <- make_speaker("spk10", n_full, 900, sd1 = 40, sd2 = 90)
spk11 <- make_speaker("spk11", 60, 110)
spk12 <- make_speaker("spk12", 40, 140)

vowel_cohort <- do.call(rbind, c(graded, list(spk09, spk10, spk11, spk12)))
rownames(vowel_cohort) <- NULL
vowel_cohort$f1 <- round(vowel_cohort$f1, 1)
vowel_cohort$f2 <- round(vowel_cohort$f2, 1)

dir.create("data", showWarnings = FALSE)
save(vowel_cohort, file = "data/vowel_cohort.rda", compress = "bzip2", version = 2)
