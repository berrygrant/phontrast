test_that("discrete JSD handles zeros correctly", {
  expect_equal(jsd(c(1, 0), c(1, 0)), 0)
  expect_equal(jsd(c(1, 0), c(0, 1)), 1)
  expect_error(jsd(c(1, -1), c(1, 1)), "finite, non-negative")
})

test_that("classical grouped metrics drop incomplete rows consistently", {
  set.seed(1)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 20),
    category = rep(rep(c("a", "b"), each = 10), 2),
    f1 = rnorm(40),
    f2 = rnorm(40)
  )
  data$f1[1] <- NA_real_

  pillai <- speaker_pillai(
    data = data,
    group_col = "speaker",
    category_col = "category",
    features = c("f1", "f2"),
    min_tokens = 10
  )
  bhatt <- estimate_bhatt(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 10
  )

  expect_equal(nrow(pillai), 2)
  expect_equal(nrow(bhatt), 2)
  expect_equal(pillai$n_tokens[pillai$group == "s1"], 19)
  expect_equal(bhatt$n_tokens[bhatt$group == "s1"], 19)
})

test_that("KDE JSD wrappers return valid global and grouped estimates", {
  set.seed(2)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )
  data$f2[1] <- NA_real_

  point <- jsd_kde_nd(data, c("f1", "f2"), "category")
  global <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    min_tokens = 20
  )
  grouped <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20
  )

  expect_true(is.finite(point))
  expect_true(point >= 0 && point <= 1)
  expect_equal(global$scope, "global")
  expect_equal(global$n_tokens, 159)
  expect_equal(nrow(grouped), 2)
  expect_true(all(is.finite(grouped$jsd_point)))
})

test_that("bootstrap JSD summaries preserve successful replicate counts", {
  set.seed(3)
  data <- data.frame(
    speaker = rep("s1", 80),
    category = rep(c("a", "b"), each = 40),
    f1 = c(rnorm(40, 0), rnorm(40, 1)),
    f2 = c(rnorm(40, 0), rnorm(40, 1))
  )

  out <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 3,
    min_tokens = 20
  )

  expect_equal(out$n_boot, 3)
  expect_equal(out$conf_level, 0.95)
  expect_true(is.finite(out$jsd_mean))
  expect_true(is.finite(out$jsd_low))
  expect_true(is.finite(out$jsd_high))
  expect_equal(out$ci_lower, out$jsd_low)
  expect_equal(out$ci_upper, out$jsd_high)
})

test_that("one-dimensional KDE metrics are supported", {
  set.seed(4)
  data <- data.frame(
    category = rep(c("a", "b"), each = 40),
    f1 = c(rnorm(40, 0), rnorm(40, 1))
  )

  jsd_val <- jsd_kde_nd(data, "f1", "category")
  overlap <- percent_overlap_kde(data, "f1", "category")

  expect_true(is.finite(jsd_val))
  expect_true(jsd_val >= 0 && jsd_val <= 1)
  expect_true(is.finite(overlap))
  expect_true(overlap >= 0 && overlap <= 1)
})

test_that("fast KDE controls support diagonal Scott bandwidths", {
  set.seed(401)
  features <- paste0("x", 1:3)
  data <- data.frame(
    category = rep(c("a", "b"), each = 60),
    x1 = c(rnorm(60, 0), rnorm(60, 0.8)),
    x2 = c(rnorm(60, 0), rnorm(60, 0.5)),
    x3 = c(rnorm(60, 0), rnorm(60, 0.3))
  )

  slow <- jsd_kde_nd(
    data,
    features,
    "category",
    bw = "scott.diag",
    eval_n = 50,
    eval_seed = 2026
  )
  fast <- jsd_kde_nd(
    data,
    features,
    "category",
    bw = "scott.diag",
    eval_n = 50,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  repeat_fast <- jsd_kde_nd(
    data,
    features,
    "category",
    bw = "scott.diag",
    eval_n = 50,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  overlap <- percent_overlap_kde(
    data,
    features,
    "category",
    bw = "scott.diag",
    eval_on = "pooled_sample",
    eval_n = 50,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  metrics <- compare_overlap_metrics(
    data = data,
    features = features,
    category_col = "category",
    bw = "scott.diag",
    eval_n = 50,
    eval_seed = 2026,
    engine = "fast_diag",
    min_tokens = 20
  )
  wrapped_fast <- estimate_jsd(
    data = data,
    features = features,
    category_col = "category",
    bw = "scott.diag",
    eval_n = 50,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  wrapped_alias <- estimate_jsd(
    data = data,
    features = features,
    category_col = "category",
    bw = "scott.diag",
    eval_n = 50,
    eval_seed = 2026,
    engine = "fast_diagonal"
  )

  expect_equal(fast, slow, tolerance = 1e-6)
  expect_equal(repeat_fast, fast)
  expect_equal(wrapped_fast$jsd_point, fast)
  expect_equal(wrapped_alias$jsd_point, fast)
  expect_true(is.finite(overlap))
  expect_true(overlap >= 0 && overlap <= 1)
  expect_equal(nrow(metrics), 1)
  expect_true(is.finite(metrics$jsd))
  expect_true(is.finite(metrics$percent_overlap))
  expect_error(
    jsd_kde_nd(data, features, "category", bw = "Hpi", engine = "fast_diag"),
    "requires `bw = \"scott.diag\"`"
  )
  expect_error(
    jsd_kde_nd(data, features, "category", eval_on = "pooled_sample"),
    "`eval_n` must be supplied"
  )
})

test_that("metric wrappers accept multiple grouping columns", {
  set.seed(402)
  data <- expand.grid(
    sex = c("F", "M"),
    style = c("casual", "read"),
    category = c("a", "b"),
    rep = seq_len(30),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  category_shift <- ifelse(data$category == "b", 1, 0)
  sex_shift <- ifelse(data$sex == "M", 0.2, 0)
  style_shift <- ifelse(data$style == "read", -0.2, 0)
  data$f1 <- rnorm(nrow(data), mean = category_shift + sex_shift)
  data$f2 <- rnorm(nrow(data), mean = category_shift + style_shift)
  group_cols <- c("sex", "style")

  jsd_out <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = group_cols,
    min_tokens = 20,
    bw = "scott.diag",
    eval_n = 30,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  summary_out <- jsd_summary(
    data = data,
    group_col = group_cols,
    category_col = "category",
    features = c("f1", "f2"),
    do_boot = FALSE,
    min_tokens = 20,
    bw = "scott.diag",
    eval_n = 30,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  overlap_out <- estimate_overlap(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = group_cols,
    min_tokens = 20,
    bw = "scott.diag",
    eval_n = 30,
    eval_seed = 2026,
    engine = "fast_diag"
  )
  pillai_out <- estimate_pillai(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = group_cols,
    min_tokens = 20
  )
  bhatt_out <- estimate_bhatt(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = group_cols,
    min_tokens = 20
  )
  metrics <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = group_cols,
    min_tokens = 20,
    bw = "scott.diag",
    eval_n = 30,
    eval_seed = 2026,
    engine = "fast_diagonal"
  )

  expected_groups <- c(
    "sex=F | style=casual",
    "sex=M | style=casual",
    "sex=F | style=read",
    "sex=M | style=read"
  )
  expect_equal(nrow(jsd_out), 4)
  expect_setequal(jsd_out$group, expected_groups)
  expect_setequal(summary_out$group, expected_groups)
  expect_setequal(overlap_out$group, expected_groups)
  expect_setequal(pillai_out$group, expected_groups)
  expect_setequal(bhatt_out$group, expected_groups)
  expect_setequal(metrics$group, expected_groups)
  expect_true(all(is.finite(jsd_out$jsd_point)))
  expect_true(all(is.finite(overlap_out$overlap)))
  expect_true(all(is.finite(metrics$jsd)))
  expect_true(all(is.finite(metrics$percent_overlap)))
})

test_that("KDE metrics ignore unused factor levels after filtering", {
  set.seed(41)
  data <- data.frame(
    category = factor(rep(c("a", "b"), each = 30), levels = c("a", "b", "unused")),
    f1 = c(rnorm(30, 0), rnorm(30, 1)),
    f2 = c(rnorm(30, 0), rnorm(30, 1))
  )

  jsd_val <- jsd_kde_nd(data, c("f1", "f2"), "category")
  overlap <- percent_overlap_kde(data, c("f1", "f2"), "category")
  metrics <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    min_tokens = 20
  )

  expect_true(is.finite(jsd_val))
  expect_true(is.finite(overlap))
  expect_equal(nrow(metrics), 1)
  expect_true(is.finite(metrics$jsd))
  expect_true(is.finite(metrics$percent_overlap))
})

test_that("global JSD wrapper forwards KDE controls", {
  set.seed(5)
  data <- data.frame(
    category = rep(c("a", "b"), each = 60),
    f1 = c(rnorm(60, 0), rnorm(60, 1)),
    f2 = c(rnorm(60, 0), rnorm(60, 1))
  )

  direct <- jsd_kde_nd(data, c("f1", "f2"), "category", bw = "Hpi.diag")
  wrapped <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    min_tokens = 10,
    bw = "Hpi.diag"
  )

  expect_equal(wrapped$jsd_point, direct)
})

test_that("grouped bootstrap JSD respects confidence level", {
  set.seed(6)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = rnorm(160),
    f2 = rnorm(160)
  )

  set.seed(7)
  ci50 <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 12,
    min_tokens = 20,
    conf_level = 0.50
  )
  set.seed(7)
  ci95 <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 12,
    min_tokens = 20,
    conf_level = 0.95
  )

  width50 <- ci50$jsd_high - ci50$jsd_low
  width95 <- ci95$jsd_high - ci95$jsd_low
  expect_equal(ci50$conf_level, rep(0.50, 2))
  expect_equal(ci95$conf_level, rep(0.95, 2))
  expect_equal(ci50$ci_lower, ci50$jsd_low)
  expect_equal(ci50$ci_upper, ci50$jsd_high)
  expect_true(all(width50 <= width95))
})

test_that("JSD bootstrap outputs include standard CI columns and confidence level", {
  set.seed(17)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )

  set.seed(18)
  global <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    do_boot = TRUE,
    n_boot = 4,
    min_tokens = 20,
    conf_level = 0.80
  )
  expect_true(all(c("conf_level", "ci_lower", "ci_upper", "jsd_low", "jsd_high") %in% names(global)))
  expect_equal(global$conf_level, 0.80)
  expect_equal(global$ci_lower, global$jsd_low)
  expect_equal(global$ci_upper, global$jsd_high)

  set.seed(18)
  grouped <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    do_boot = TRUE,
    n_boot = 4,
    min_tokens = 20,
    conf_level = 0.80
  )
  expect_true(all(c("conf_level", "ci_lower", "ci_upper", "jsd_low", "jsd_high") %in% names(grouped)))
  expect_equal(grouped$conf_level, rep(0.80, 2))
  expect_equal(grouped$ci_lower, grouped$jsd_low)
  expect_equal(grouped$ci_upper, grouped$jsd_high)

  no_boot <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20
  )
  expect_equal(no_boot$n_boot, rep(0L, 2))
  expect_equal(no_boot$conf_level, rep(0.95, 2))
  expect_true(all(is.na(no_boot$ci_lower)))
  expect_true(all(is.na(no_boot$ci_upper)))
})

test_that("lower-level JSD summaries carry standard CI columns", {
  set.seed(19)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )

  set.seed(20)
  boot <- boot_jsd(
    data = data,
    group_col = "speaker",
    category_col = "category",
    features = c("f1", "f2"),
    n_boot = 4,
    min_tokens = 20,
    conf_level = 0.90
  )
  expect_equal(boot$conf_level, rep(0.90, 2))
  expect_equal(boot$ci_lower, boot$jsd_low)
  expect_equal(boot$ci_upper, boot$jsd_high)

  no_boot_summary <- jsd_summary(
    data = data,
    group_col = "speaker",
    category_col = "category",
    features = c("f1", "f2"),
    do_boot = FALSE,
    min_tokens = 20
  )
  expect_equal(no_boot_summary$n_boot, rep(0L, 2))
  expect_equal(no_boot_summary$conf_level, rep(0.95, 2))
  expect_true(all(is.na(no_boot_summary$ci_lower)))
  expect_true(all(is.na(no_boot_summary$ci_upper)))

  set.seed(21)
  global <- global_boot_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    n_boot = 4,
    min_tokens = 20,
    conf_level = 0.90
  )
  expect_equal(global$conf_level, 0.90)
  expect_equal(global$ci_lower, global$jsd_low)
  expect_equal(global$ci_upper, global$jsd_high)
})

test_that("grouped metrics return shaped empty outputs when all groups are filtered", {
  set.seed(8)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 20),
    category = rep(rep(c("a", "b"), each = 10), 2),
    f1 = rnorm(40),
    f2 = rnorm(40)
  )

  pillai <- estimate_pillai(data, c("f1", "f2"), "category", "speaker", min_tokens = 100)
  bhatt <- estimate_bhatt(data, c("f1", "f2"), "category", "speaker", min_tokens = 100)
  overlap <- estimate_overlap(data, c("f1", "f2"), "category", "speaker", min_tokens = 100)

  expect_equal(nrow(pillai), 0)
  expect_equal(names(pillai), c("scope", "group", "n_tokens", "pillai", "p_value"))
  expect_equal(nrow(bhatt), 0)
  expect_equal(names(bhatt), c("scope", "group", "n_tokens", "bhatt_dist", "bhatt_affinity"))
  expect_equal(nrow(overlap), 0)
  expect_equal(names(overlap), c("scope", "group", "n_tokens", "overlap"))
})

test_that("prepare_jsd_beta rejects invalid JSD values", {
  ok <- prepare_jsd_beta(data.frame(jsd_mean = c(0, 0.5, 1, NA)), eps = 0.01)
  expect_equal(ok$jsd_beta, c(0.01, 0.5, 0.99, NA))
  expect_error(
    prepare_jsd_beta(data.frame(jsd_mean = c(-0.1, 0.2, Inf))),
    "values in \\[0, 1\\]"
  )
})

test_that("small KDE samples fail with a clear data-size message", {
  set.seed(9)
  data <- data.frame(
    category = rep(c("a", "b"), each = 2),
    f1 = rnorm(4),
    f2 = rnorm(4)
  )

  expect_error(
    jsd_kde_nd(data, c("f1", "f2"), "category"),
    "at least 3 finite observations"
  )
})

test_that("compare_overlap_metrics returns wide and long comparisons", {
  set.seed(10)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )

  wide <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20
  )
  long <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20,
    output = "long"
  )

  expect_equal(nrow(wide), 2)
  expect_true(all(c(
    "pillai", "bhatt_dist", "bhatt_affinity", "jsd", "js_distance",
    "mahalanobis_dist", "percent_overlap"
  ) %in% names(wide)))
  expect_equal(nrow(long), 14)
  expect_true(all(c(
    "metric", "estimate", "orientation", "separation_value", "separation_rank"
  ) %in% names(long)))
})

test_that("metric plotting accepts wide and long comparison output", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(201)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 60),
    category = rep(rep(c("a", "b"), each = 30), 2),
    f1 = c(rnorm(30, 0), rnorm(30, 1), rnorm(30, 0.2), rnorm(30, 1.2)),
    f2 = c(rnorm(30, 0), rnorm(30, 1), rnorm(30, 0.2), rnorm(30, 1.2))
  )

  wide <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20,
    output = "wide"
  )
  long <- compare_overlap_metrics(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20,
    output = "long"
  )

  expect_s3_class(plot_overlap_metrics(wide), "ggplot")
  expect_s3_class(plot_overlap_metrics(long, value = "estimate", facet = FALSE), "ggplot")
})

test_that("category space plotting supports one and two dimensions", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(202)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = factor(rep(rep(c("a", "b"), each = 40), 2), levels = c("a", "b", "unused")),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )

  expect_s3_class(plot_category_space(data, "f1", "category"), "ggplot")
  expect_s3_class(
    plot_category_space(
      data,
      c("f1", "f2"),
      "category",
      group_col = "speaker",
      ellipses = FALSE
    ),
    "ggplot"
  )
  expect_error(
    plot_category_space(data, c("f1", "f2", "f3"), "category"),
    "one or two column names"
  )
})

test_that("PCA category plotting supports multidimensional feature sets", {
  testthat::skip_if_not_installed("ggplot2")
  set.seed(203)
  features <- paste0("x", 1:5)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 50),
    category = rep(c("a", "b"), each = 50),
    matrix(rnorm(100 * length(features)), ncol = length(features))
  )
  names(data)[-(1:2)] <- features
  data[data$category == "b", features[1:2]] <- data[data$category == "b", features[1:2]] + 0.75

  p <- plot_category_pca(
    data = data,
    features = features,
    category_col = "category",
    group_col = "speaker"
  )

  expect_s3_class(p, "ggplot")
  expect_s3_class(attr(p, "pca"), "prcomp")
  expect_true(all(c("component", "variance_explained") %in% names(attr(p, "variance_explained"))))
  expect_error(
    plot_category_pca(data, features = "x1", category_col = "category"),
    "at least two"
  )
})

test_that("compare_overlap_metrics diagnoses grouped contrasts with no estimable groups", {
  set.seed(21)
  data <- data.frame(
    speaker = rep(paste0("s", 1:5), each = 4),
    category = factor(rep(c("I", "I", "i", "i"), 5), levels = c("i", "I", "E")),
    f1 = rnorm(20),
    f2 = rnorm(20)
  )

  expect_warning(
    metrics <- compare_overlap_metrics(
      data = data,
      features = c("f1", "f2"),
      category_col = "category",
      group_col = "speaker"
    ),
    "returned no grouped rows"
  )
  expect_equal(nrow(metrics), 0)
})

test_that("compare_overlap_metrics can bootstrap all reported metrics", {
  set.seed(22)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 60),
    category = rep(rep(c("a", "b"), each = 30), 2),
    f1 = c(rnorm(30, 0), rnorm(30, 1), rnorm(30, 0.2), rnorm(30, 1.2))
  )

  set.seed(23)
  wide <- compare_overlap_metrics(
    data = data,
    features = "f1",
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20,
    do_boot = TRUE,
    n_boot = 2,
    conf_level = 0.80,
    progress = FALSE,
    output = "wide"
  )

  expect_equal(nrow(wide), 2)
  expect_equal(wide$n_boot, rep(2, 2))
  expect_equal(wide$conf_level, rep(0.80, 2))
  expect_true(all(c(
    "jsd_mean", "jsd_sd", "jsd_ci_lower", "jsd_ci_upper", "jsd_n_boot",
    "percent_overlap_mean", "percent_overlap_ci_lower", "percent_overlap_ci_upper",
    "mahalanobis_dist_mean", "bhatt_affinity_mean"
  ) %in% names(wide)))
  expect_true(all(wide$jsd_n_boot <= 2))
  expect_true(all(is.finite(wide$jsd_ci_lower)))
  expect_true(all(is.finite(wide$jsd_ci_upper)))

  set.seed(23)
  long <- compare_overlap_metrics(
    data = data,
    features = "f1",
    category_col = "category",
    group_col = "speaker",
    min_tokens = 20,
    do_boot = TRUE,
    n_boot = 2,
    conf_level = 0.80,
    progress = FALSE,
    output = "long"
  )

  expect_equal(nrow(long), 14)
  expect_true(all(c(
    "boot_mean", "boot_sd", "ci_lower", "ci_upper", "n_boot", "conf_level"
  ) %in% names(long)))
  expect_equal(unique(long$conf_level), 0.80)
  expect_true(all(long$n_boot <= 2))
})

test_that("compare_overlap_metrics reports progress while bootstrapping", {
  set.seed(24)
  data <- data.frame(
    category = rep(c("a", "b"), each = 30),
    f1 = c(rnorm(30, 0), rnorm(30, 1))
  )

  messages <- character()
  withCallingHandlers(
    compare_overlap_metrics(
      data = data,
      features = "f1",
      category_col = "category",
      min_tokens = 20,
      do_boot = TRUE,
      n_boot = 2,
      progress = TRUE
    ),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("Bootstrapping overlap metrics", messages)))
  expect_true(any(grepl("This may take time", messages)))
  expect_true(any(grepl("bootstrap replicate 1 / 2", messages)))
  expect_true(any(grepl("bootstrap replicate 2 / 2", messages)))
})

test_that("hierarchical bootstrap preserves repeated group draws", {
  set.seed(11)
  data <- data.frame(
    speaker = rep(c("s1", "s2"), each = 80),
    category = rep(rep(c("a", "b"), each = 40), 2),
    f1 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2)),
    f2 = c(rnorm(40, 0), rnorm(40, 1), rnorm(40, 0.2), rnorm(40, 1.2))
  )
  rows_seen <- integer()
  fit_fun <- function(formula, data, ...) {
    rows_seen <<- c(rows_seen, nrow(data))
    stats::lm(formula = formula, data = data, ...)
  }

  set.seed(12)
  out <- hier_boot_jsd_model(
    data = data,
    group_col = "speaker",
    category_col = "category",
    features = c("f1", "f2"),
    formula = jsd_beta ~ 1,
    fit_fun = fit_fun,
    n_outer = 4,
    min_tokens = 20,
    progress = FALSE
  )

  expect_true(nrow(out) > 0)
  expect_true(all(rows_seen == 2))
})

test_that("extract_mfcc validates arguments before optional audio dependencies", {
  expect_error(
    extract_mfcc(data.frame(file = "missing.wav"), "file", numcep = 0),
    "positive integer"
  )
  expect_error(
    extract_mfcc(data.frame(file = "missing.wav"), "file", fs = -1),
    "positive finite"
  )
  expect_error(
    extract_mfcc(data.frame(file = "missing.wav", start = "0", end = 1), "file",
                 start_col = "start", end_col = "end"),
    "numeric column"
  )
})

test_that("metric entry points validate their inputs", {
  data <- data.frame(
    category = rep(c("a", "b"), each = 20),
    f1 = rnorm(40),
    f2 = rnorm(40)
  )

  expect_error(
    estimate_jsd(as.matrix(data[, c("f1", "f2")]), features = "f1",
                 category_col = "category"),
    "must be a data frame"
  )
  expect_error(
    estimate_jsd(data, features = c("f1", "f2"),
                 category_col = c("category", "f1")),
    "single column name"
  )
  expect_error(
    compare_overlap_metrics(data, features = c("f1", "category"),
                            category_col = "category"),
    "must not overlap"
  )
})

test_that("Jensen-Shannon distance stays finite and non-negative near zero", {
  set.seed(778)
  data <- data.frame(
    category = rep(c("a", "b"), each = 30),
    f1 = c(rnorm(30, 0), rnorm(30, 0.05)),
    f2 = c(rnorm(30, 0), rnorm(30, 0.05))
  )

  out <- estimate_jsd(
    data = data,
    features = c("f1", "f2"),
    category_col = "category",
    est_distance = TRUE,
    min_tokens = 20
  )

  expect_true(is.finite(out$jsd_point))
  expect_true(out$jsd_point >= 0)
})

test_that("grouped metrics keep failed groups as NA and warn consistently", {
  set.seed(777)
  # The "bad" speaker passes the min_tokens/two-category filter (4 tokens, both
  # categories present) but has too few tokens in category "b" for KDE/covariance
  # estimation, so the per-group metric computation fails.
  data <- data.frame(
    speaker = c(rep("good", 40), rep("bad", 4)),
    category = c(rep(c("a", "b"), each = 20), "a", "a", "a", "b"),
    f1 = rnorm(44),
    f2 = rnorm(44)
  )

  expect_warning(
    bhatt <- estimate_bhatt(
      data = data,
      features = c("f1", "f2"),
      category_col = "category",
      group_col = "speaker",
      min_tokens = 4
    ),
    "could not be estimated"
  )
  expect_warning(
    jsd <- estimate_jsd(
      data = data,
      features = c("f1", "f2"),
      category_col = "category",
      group_col = "speaker",
      min_tokens = 4
    ),
    "could not be estimated"
  )

  # Both wrappers keep the same set of groups (the failed one is retained as NA,
  # not silently dropped), so row counts agree across metrics.
  expect_setequal(bhatt$group, c("good", "bad"))
  expect_setequal(jsd$group, bhatt$group)
  expect_true(is.na(bhatt$bhatt_dist[bhatt$group == "bad"]))
  expect_true(is.finite(bhatt$bhatt_dist[bhatt$group == "good"]))
  expect_true(is.na(jsd$jsd_point[jsd$group == "bad"]))
})
