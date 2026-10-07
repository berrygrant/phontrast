# Regression test: rank_contrasts() reproduces the published calibration of the
# ranking protocol (Berry, under review, Sec. VII A) on the 45 Peterson--Barney
# F1 x F2 vowel pairs. The fixtures and their provenance are described in
# fixtures/README.md. The paper computed sqrt(JSD) with the legacy
# (self-normalized) kernel estimator at the plug-in bandwidth, and the
# bandwidth check at half and twice one diagonal Scott bandwidth selected on
# the pooled pair; phontrast 2.5.0's Monte-Carlo estimator put 20 pairs at the
# ceiling instead of 23 and set no pair aside.

pb52_tokens <- utils::read.csv(
  test_path("fixtures", "pb52_f1f2.csv"),
  stringsAsFactors = FALSE, encoding = "UTF-8"
)
pb52_paper <- utils::read.csv(
  test_path("fixtures", "pb52_paper_per_pair.csv"),
  stringsAsFactors = FALSE, encoding = "UTF-8"
)
pb52_paper$group <- paste(pb52_paper$label_a, pb52_paper$label_b, sep = " | ")

# One "speaker" per vowel pair, so each group is one measurement of one
# contrast, as in the paper's Peterson--Barney walk-through.
pb52_stacked <- do.call(rbind, lapply(seq_len(nrow(pb52_paper)), function(i) {
  d <- pb52_tokens[pb52_tokens$vowel %in% c(pb52_paper$label_a[i], pb52_paper$label_b[i]), ]
  d$pair <- pb52_paper$group[i]
  d
}))

pb52_ranking <- rank_contrasts(pb52_stacked, c("f1", "f2"), "vowel", "pair")
pb52_merged <- merge(as.data.frame(pb52_ranking), pb52_paper, by = "group")

test_that("the PB52 fixtures are the paper's 45 pairs of 152 tokens", {
  expect_equal(nrow(pb52_tokens), 1520L)
  expect_true(all(table(pb52_tokens$vowel) == 152L))
  expect_equal(nrow(pb52_paper), 45L)
  expect_equal(nrow(pb52_merged), 45L)
})

test_that("rank_contrasts() reproduces the paper's PB52 outcome", {
  r <- pb52_ranking
  p <- attr(r, "protocol")
  expect_identical(p$estimator$method, "legacy")
  expect_identical(p$estimator$bracket_bw, "scott.pooled")
  expect_equal(nrow(r), 45L)
  expect_equal(sum(r$at_ceiling), 23L)
  expect_equal(p$k, 22L)
  expect_equal(sum(!is.na(r$pr_jsd)), 22L)
  expect_equal(sum(r$flag %in% TRUE), 0L)
  expect_identical(r$group[r$set_aside %in% TRUE], "3' | u")
  # ɝ-u is set aside by a rank shift of 6/22 against the 0.25 margin
  expect_equal(abs(r$bw_shift[r$group == "3' | u"]), 6 / 22)
  expect_true(all(r$rank_licensed))
  expect_true(all(r$flag_licensed))
})

test_that("rank_contrasts() reproduces the paper's PB52 per-pair values", {
  m <- pb52_merged
  lg <- function(x) as.character(x) %in% c("TRUE", "True")
  # step 1: sqrt(JSD), shared mass, and Pillai on the same tokens
  expect_equal(m$sqrt_jsd, m$step1_sqrt_jsd, tolerance = 1e-10)
  expect_equal(m$shared_mass, m$step1_kde_overlap, tolerance = 1e-10)
  expect_equal(m$pillai, m$step1_pillai, tolerance = 1e-10)
  # ɑ-ɔ (internal label A | O), where the Monte-Carlo estimator read 0.730
  expect_equal(round(m$sqrt_jsd[m$group == "A | O"], 3), 0.627)
  # the ceiling
  expect_identical(m$at_ceiling, lg(m$step5_saturated))
  # steps 2 and 3, over the 22 pairs below the ceiling
  b <- !m$at_ceiling
  expect_equal(m$pr_jsd[b], m$step2_pct_sqrt_jsd_Bprime[b], tolerance = 1e-12)
  expect_equal(m$pr_pillai[b], m$step2_pct_pillai_Bprime[b], tolerance = 1e-12)
  expect_identical(m$flag[b], lg(m$step2_flag_Bprime[b]))
  # the bandwidth check
  expect_equal(m$sqrt_jsd_half, m$step4_sqrt_jsd_x0p5, tolerance = 1e-10)
  expect_equal(m$sqrt_jsd_double, m$step4_sqrt_jsd_x2, tolerance = 1e-10)
  expect_equal(abs(m$bw_shift[b]), m$step4_rank_shift_Bprime[b], tolerance = 1e-12)
  expect_identical(m$set_aside[b], lg(m$step4_set_aside_Bprime[b]))
})
