# Tests for rank_contrasts(), the one-call Sec. VII.A protocol, and for
# protocol_floors().

# A cohort of eight speakers with a graded F1/F2 mean gap, one planted
# disagreement (equal means, so Pillai ~ 0, but a bimodal second category, so
# sqrt(JSD) is substantial), and one speaker at the ceiling.
cohort_fixture <- function(n = 100, seed = 11, gaps = seq(25, 200, length.out = 8),
                           plant = TRUE, ceiling = TRUE) {
  set.seed(seed)
  make <- function(id, gap, sd1 = 55, sd2 = 120) {
    data.frame(
      speaker = id,
      vowel = rep(c("ih", "eh"), each = n),
      f1 = c(rnorm(n, 500, sd1), rnorm(n, 500 + gap, sd1)),
      f2 = c(rnorm(n, 1900, sd2), rnorm(n, 1900 - gap, sd2)),
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, Map(make, sprintf("s%02d", seq_along(gaps)), gaps))
  if (plant) {
    side <- sample(c(-1, 1), n, replace = TRUE)
    out <- rbind(out, data.frame(
      speaker = "plant",
      vowel = rep(c("ih", "eh"), each = n),
      f1 = c(rnorm(n, 500, 30), 500 + 95 * side + rnorm(n, 0, 30)),
      f2 = c(rnorm(n, 1900, 70), 1900 - 220 * side + rnorm(n, 0, 70)),
      stringsAsFactors = FALSE
    ))
  }
  if (ceiling) {
    out <- rbind(out, make("ceil", 900, sd1 = 40, sd2 = 90))
  }
  out
}

# The default-settings ranking of the full cohort, computed once and shared.
cohort_main <- cohort_fixture()
ranking_main <- rank_contrasts(cohort_main, c("f1", "f2"), "vowel", "speaker")

# ---- protocol_floors() -------------------------------------------------------

test_that("protocol_floors() encodes the Sec. V.D / VII.A floors", {
  f2 <- protocol_floors(2)
  expect_equal(f2$rank_floor, 50)
  expect_equal(f2$flag_floor, 100)
  expect_false(f2$ordering_only)
  expect_equal(protocol_floors(1)$rank_floor, 50)
  f4 <- protocol_floors(4)
  expect_equal(f4$rank_floor, 200)
  expect_equal(f4$flag_floor, 200)
  # untabulated dimensionalities take the next higher calibrated floors
  expect_equal(protocol_floors(3), modifyList(f4, list(d = 3L)))
  f8 <- protocol_floors(8)
  expect_equal(f8$rank_floor, 500)
  expect_identical(f8$flag_floor, Inf)
  expect_equal(protocol_floors(5)$rank_floor, 500)
  f13 <- protocol_floors(13)
  expect_identical(f13$rank_floor, Inf)
  expect_true(f13$ordering_only)
  expect_error(protocol_floors(0), "positive integer")
  expect_error(protocol_floors(c(2, 4)), "positive integer")
})

# ---- rank_contrasts(): structure and the three steps -------------------------

test_that("rank_contrasts() returns the protocol table, sorted by sqrt(JSD)", {
  d <- cohort_main
  r <- ranking_main
  expect_s3_class(r, "phontrast_ranking")
  expect_s3_class(r, "tbl_df")
  expect_equal(nrow(r), 10L)
  expect_named(r, c(
    "group", "n_tokens", "n_min", "sqrt_jsd", "pillai", "shared_mass",
    "at_ceiling", "pr_jsd", "pr_pillai", "rank_diff", "flag", "rank_licensed",
    "flag_licensed", "rank_basis", "sqrt_jsd_half", "sqrt_jsd_double",
    "bw_shift", "sign_change", "set_aside"
  ))
  expect_equal(r$sqrt_jsd, sort(r$sqrt_jsd, decreasing = TRUE))
  expect_true(all(r$n_min == 100L))
  expect_true(all(r$n_tokens == 200L))
  expect_true(all(r$sqrt_jsd >= 0 & r$sqrt_jsd <= 1))
  expect_true(all(r$shared_mass >= 0 & r$shared_mass <= 1))
  expect_true(all(r$pillai >= 0 & r$pillai <= 1))
  p <- attr(r, "protocol")
  expect_equal(p$d, 2L)
  expect_equal(p$margin, 0.25)
  expect_equal(p$estimator$bw, "Hpi")
})

test_that("step 1 matches the package's own estimators on the same tokens", {
  d <- cohort_fixture(n = 60, plant = FALSE, ceiling = FALSE)
  r <- suppressWarnings(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                                       bw_check = FALSE))
  s03 <- d[d$speaker == "s03", ]
  expect_equal(r$sqrt_jsd[r$group == "s03"],
               sqrt(jsd_kde_nd(s03, c("f1", "f2"), "vowel")))
  expect_equal(r$shared_mass[r$group == "s03"],
               percent_overlap_kde(s03, c("f1", "f2"), "vowel"))
  expect_equal(r$pillai[r$group == "s03"],
               pillai_overlap(s03, c("f1", "f2"), "vowel")$pillai)
})

test_that("step 2 ranks on the percentile scale and recovers a graded ordering", {
  d <- cohort_main
  r <- ranking_main
  ranked <- r[!r$at_ceiling, ]
  k <- nrow(ranked)
  expect_equal(sort(ranked$pr_jsd), (seq_len(k) - 0.5) / k)
  expect_equal(sort(ranked$pr_pillai), (seq_len(k) - 0.5) / k)
  expect_equal(ranked$rank_diff, ranked$pr_jsd - ranked$pr_pillai)
  graded <- ranked[grepl("^s", ranked$group), ]
  gap <- as.integer(sub("s", "", graded$group))
  expect_gt(stats::cor(graded$sqrt_jsd, gap, method = "spearman"), 0.9)
  expect_gt(stats::cor(graded$pillai, gap, method = "spearman"), 0.9)
})

test_that("the ceiling speaker is set apart but keeps Pillai and shared mass", {
  d <- cohort_main
  r <- ranking_main
  ceil <- r[r$group == "ceil", ]
  expect_true(ceil$at_ceiling)
  expect_gte(ceil$sqrt_jsd, 0.99)
  expect_true(is.na(ceil$pr_jsd))
  expect_true(is.na(ceil$pr_pillai))
  expect_true(is.na(ceil$flag))
  expect_false(is.na(ceil$pillai))
  expect_false(is.na(ceil$shared_mass))
  expect_equal(attr(r, "protocol")$k, 9L)
  expect_equal(attr(r, "protocol")$n_ceiling, 1L)
  # a stricter ceiling pulls more speakers out of the ranking
  r2 <- suppressWarnings(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                                        ceiling = 0.5, bw_check = FALSE))
  expect_gt(sum(r2$at_ceiling), 1L)
  expect_equal(sum(!is.na(r2$pr_jsd)), sum(!r2$at_ceiling))
})

test_that("step 3 flags the planted Pillai / sqrt(JSD) disagreement", {
  d <- cohort_main
  r <- ranking_main
  plant <- r[r$group == "plant", ]
  expect_lt(plant$pillai, 0.1)
  expect_gt(plant$sqrt_jsd, 0.4)
  expect_gte(plant$rank_diff, 0.25)
  expect_true(plant$flag)
  expect_true(plant$flag_licensed)
  expect_true(plant$rank_licensed)
  expect_identical(plant$rank_basis, "sqrt_jsd")
  # the Gaussian speakers agree
  expect_false(any(r$flag[grepl("^s", r$group)]))
  # the margin is a dial
  r_loose <- rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                            margin = 0.9, bw_check = FALSE)
  expect_false(any(r_loose$flag, na.rm = TRUE))
})

# ---- licensing ---------------------------------------------------------------

test_that("floors are applied per speaker from the smaller category's count", {
  d <- cohort_fixture(n = 30, plant = FALSE, ceiling = FALSE)
  # at d = 2: rank from 50, flag from 100 -> neither licensed at 30 tokens
  r <- suppressWarnings(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                                       bw_check = FALSE))
  expect_true(all(!r$rank_licensed))
  expect_true(all(r$rank_basis == "pillai"))
  expect_true(all(!r$flag_licensed))
  expect_true(all(is.na(r$flag)))
  expect_false(any(is.na(r$rank_diff)))

  # 60 tokens: ranking licensed, flag still not readable
  d60 <- cohort_fixture(n = 60, plant = FALSE, ceiling = FALSE)
  r60 <- suppressWarnings(rank_contrasts(d60, c("f1", "f2"), "vowel", "speaker",
                                         bw_check = FALSE))
  expect_true(all(r60$rank_licensed))
  expect_true(all(!r60$flag_licensed))

  # unequal counts: the smaller category governs
  d_uneq <- cohort_fixture(n = 60, plant = FALSE, ceiling = FALSE)
  drop <- which(d_uneq$speaker == "s01" & d_uneq$vowel == "eh")[1:20]
  d_uneq <- d_uneq[-drop, ]
  r_uneq <- suppressWarnings(rank_contrasts(d_uneq, c("f1", "f2"), "vowel",
                                            "speaker", bw_check = FALSE))
  expect_equal(r_uneq$n_min[r_uneq$group == "s01"], 40L)
  expect_false(r_uneq$rank_licensed[r_uneq$group == "s01"])
  expect_true(all(r_uneq$rank_licensed[r_uneq$group != "s01"]))
})

test_that("three features take the d = 4 floors and the d >= 2n rule applies", {
  d <- cohort_fixture(n = 30, plant = FALSE, ceiling = FALSE)
  d$f3 <- d$f1 * 0.3 + rnorm(nrow(d), 0, 40)
  r <- suppressWarnings(rank_contrasts(d, c("f1", "f2", "f3"), "vowel", "speaker",
                                       bw_check = FALSE))
  expect_equal(attr(r, "protocol")$floors$rank_floor, 200)
  expect_true(all(!r$rank_licensed))
  expect_true(all(is.na(r$flag)))
  # the flag is undefined when d >= 2 * n_min even if the floors were met
  expect_false(phontrast:::.flag_defined(d = 4L, n_min = 2L))
  expect_true(phontrast:::.flag_defined(d = 2L, n_min = 100L))
})

# ---- bandwidth check -----------------------------------------------------------

test_that("the bandwidth check re-ranks at half and twice the Scott bandwidth", {
  d <- cohort_main
  r <- ranking_main
  ok <- !r$at_ceiling
  expect_true(all(is.finite(r$sqrt_jsd_half[ok])))
  expect_true(all(is.finite(r$sqrt_jsd_double[ok])))
  # the two re-estimates differ for every ranked speaker (their direction at
  # two dimensions depends on the sample; see test-protocol-helpers.R)
  expect_true(all(r$sqrt_jsd_half[ok] != r$sqrt_jsd_double[ok]))
  s05 <- d[d$speaker == "s05", ]
  expect_equal(
    r$sqrt_jsd_half[r$group == "s05"],
    sqrt(jsd_kde_nd(s05, c("f1", "f2"), "vowel", bw = "scott.diag", bw_scale = 0.5))
  )
  expect_true(is.logical(r$set_aside))
  expect_true(is.logical(r$sign_change))
  expect_true(all(abs(r$bw_shift[ok]) <= 1))
  expect_true(is.na(r$bw_shift[!ok]))
  # a graded Gaussian cohort is stable under smoothing
  expect_false(any(r$set_aside[grepl("^s", r$group)]))

  r_off <- rank_contrasts(d, c("f1", "f2"), "vowel", "speaker", bw_check = FALSE)
  expect_false(any(c("sqrt_jsd_half", "sqrt_jsd_double", "bw_shift",
                     "sign_change", "set_aside") %in% names(r_off)))
})

# ---- estimator settings, dropping, errors, print -------------------------------

test_that("estimator settings are validated and honoured", {
  d <- cohort_fixture(n = 60, plant = FALSE, ceiling = FALSE)
  est <- modifyList(recommended_estimator(2), list(bw = "scott.diag", engine = "fast_diag"))
  r <- suppressWarnings(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                                       estimator = est, bw_check = FALSE))
  s02 <- d[d$speaker == "s02", ]
  expect_equal(
    r$sqrt_jsd[r$group == "s02"],
    sqrt(jsd_kde_nd(s02, c("f1", "f2"), "vowel", bw = "scott.diag", engine = "fast_diag"))
  )
  expect_equal(attr(r, "protocol")$estimator$engine, "fast_diag")
  expect_error(
    rank_contrasts(d, c("f1", "f2"), "vowel", "speaker", estimator = "Hpi"),
    "estimator"
  )
  expect_error(
    rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                   estimator = list(bw = "nope", engine = "ks", eval_n = NULL, loo = TRUE)),
    "estimator\\$bw"
  )
})

test_that("speakers with too few tokens are left out with a message", {
  d <- cohort_fixture(n = 60, plant = FALSE, ceiling = FALSE)
  d <- rbind(d, data.frame(speaker = "tiny", vowel = rep(c("ih", "eh"), each = 4),
                           f1 = rnorm(8, 500, 50), f2 = rnorm(8, 1900, 100),
                           stringsAsFactors = FALSE))
  d <- rbind(d, data.frame(speaker = "onecat", vowel = "ih",
                           f1 = rnorm(30, 500, 50), f2 = rnorm(30, 1900, 100),
                           stringsAsFactors = FALSE))
  expect_message(
    r <- suppressWarnings(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker",
                                         bw_check = FALSE)),
    "not measured: onecat, tiny|not measured: tiny, onecat"
  )
  expect_false(any(c("tiny", "onecat") %in% r$group))
  expect_equal(attr(r, "protocol")$n_dropped, 2L)
})

test_that("rank_contrasts() warns on small sets and errors when nothing can be ranked", {
  d <- cohort_fixture(n = 60, gaps = c(50, 100, 150), plant = FALSE, ceiling = FALSE)
  expect_warning(
    rank_contrasts(d, c("f1", "f2"), "vowel", "speaker", bw_check = FALSE),
    "smallest useful set"
  )
  one <- d[d$speaker == "s01", ]
  expect_error(
    rank_contrasts(one, c("f1", "f2"), "vowel", "speaker", bw_check = FALSE),
    "at least two"
  )
  expect_error(rank_contrasts(d, c("f1", "f2"), "vowel"), "group_col")
  expect_error(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker", margin = 0), "margin")
  expect_error(rank_contrasts(d, c("f1", "f2"), "vowel", "speaker", ceiling = 1.5), "ceiling")
})

test_that("rank_contrasts() works on one feature and with multiple grouping columns", {
  d <- cohort_fixture(n = 60, plant = FALSE, ceiling = FALSE)
  r1 <- suppressWarnings(rank_contrasts(d, "f1", "vowel", "speaker", bw_check = FALSE))
  expect_equal(nrow(r1), 8L)
  expect_equal(attr(r1, "protocol")$floors$rank_floor, 50)
  d$style <- "read"
  r2 <- suppressWarnings(rank_contrasts(d, c("f1", "f2"), "vowel",
                                        c("speaker", "style"), bw_check = FALSE))
  expect_true(all(grepl("speaker=.* \\| style=read", r2$group)))
})

test_that("print.phontrast_ranking() summarises the protocol", {
  d <- cohort_main
  r <- ranking_main
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_match(out, "phontrast ranking: 10 speakers, 2 features")
  expect_match(out, "1 at ceiling")
  expect_match(out, "Flagged .* -> plant")
  expect_match(out, "Bandwidth check")
  expect_match(out, "sqrt_jsd")
})
