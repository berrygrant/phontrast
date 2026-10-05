#' Rank speakers' contrasts by Jensen-Shannon distance and check them against Pillai
#'
#' One-call implementation of the measurement protocol recommended by the
#' simulation study behind phontrast (Berry, under review, Sec. VII.A). For
#' each speaker (or other unit named by \code{group_col}) it
#' \enumerate{
#'   \item computes Jensen-Shannon distance (\eqn{\sqrt{JSD}}) and the Pillai
#'     trace on the same tokens, from one shared kernel density estimate, and
#'     reports the estimated shared probability mass beside them;
#'   \item ranks the speakers by \eqn{\sqrt{JSD}} on the percentile-rank scale
#'     of \code{percentile_rank()};
#'   \item flags speakers whose Pillai percentile rank differs from their
#'     \eqn{\sqrt{JSD}} percentile rank by \code{margin} (0.25 of the ordering)
#'     or more, for inspection with \code{inspect_contrast()} or
#'     \code{plot_contrast()}.
#' }
#' The conditions the study attaches to those steps are applied and reported
#' rather than left to the user: measurements at the \eqn{\sqrt{JSD}} ceiling
#' are set apart, sample-size floors decide whether a rank or a flag may be
#' read, and a bandwidth check marks measurements whose rank depends on the
#' smoothing.
#'
#' @section Ceiling:
#' A measurement with \eqn{\sqrt{JSD} \ge} \code{ceiling} (default 0.99) is at
#' the ceiling: the kernel estimate can no longer separate degrees of
#' separation. Its shared mass and Pillai value are reported, but it is
#' excluded from both percentile rankings (\code{at_ceiling = TRUE},
#' \code{pr_jsd} and \code{pr_pillai} are \code{NA}), and the ranks of the
#' remaining \eqn{k} speakers are computed over those \eqn{k} alone.
#'
#' @section Sample-size floors:
#' Floors are read from the smaller category's token count per speaker
#' (\code{n_min}) and from the number of features \code{d}
#' (see \code{protocol_floors()}). Ranking by \eqn{\sqrt{JSD}} is licensed from
#' 50 tokens per category at \code{d = 2}, 200 at \code{d = 3} or \code{4},
#' and 500 at \code{d = 5} to \code{8}; above eight dimensions the study gives
#' ordering evidence only, so no rank is licensed. Below the floor
#' \code{rank_basis} is \code{"pillai"}: order that speaker by Pillai. The
#' flag has floors of its own: 100 tokens per category at \code{d = 2} and
#' 200 at \code{d = 3} or \code{4}; from eight dimensions up the agreement test
#' returns agreement whatever the data contain, and it is undefined when
#' \code{d >= 2 * n_min}. Where the flag may not be read, \code{flag} is
#' \code{NA} and \code{rank_diff} is still reported. The floors are where the
#' simulation recovers the average ordering of its separation levels; they do
#' not guarantee an accurate ranking of an individual speaker. Eight speakers
#' is the smallest useful set (the calibration used sets of 32); fewer draws a
#' warning.
#'
#' @section Bandwidth check:
#' With \code{bw_check = TRUE} (the default) \eqn{\sqrt{JSD}} is recomputed at
#' half and at twice the diagonal Scott bandwidth (\code{bw = "scott.diag"},
#' \code{bw_scale = 0.5} and \code{2}), at two dimensions as well, where the
#' reported estimate uses the plug-in rule. The two re-estimates are ranked
#' over the same speakers and \code{bw_shift} is the percentile-rank change
#' between them. A measurement is set aside (\code{set_aside = TRUE}) when
#' \code{abs(bw_shift) >= margin}, or when a flagged Pillai--\eqn{\sqrt{JSD}}
#' rank difference changes sign between the halved and the doubled bandwidth
#' (\code{sign_change}). The check tests sensitivity to smoothing, not the
#' estimate's accuracy.
#'
#' @param data Data frame with the speaker, category, and feature columns.
#' @param features Character vector of numeric feature columns.
#' @param category_col String; column giving the two categories of the
#'   contrast (for example the two vowels).
#' @param group_col Character vector of one or more columns identifying the
#'   speakers (or other units) to rank. Required: one speaker measured once on
#'   one contrast gives nothing to rank.
#' @param margin Inspection margin on the percentile-rank scale: the flag is
#'   raised when \code{abs(rank_diff) >= margin} and a measurement is set aside
#'   by the bandwidth check when its \eqn{\sqrt{JSD}} rank moves by
#'   \code{margin} or more. Default 0.25 of the ordering, the value calibrated
#'   in the study.
#' @param ceiling \eqn{\sqrt{JSD}} at or above which a measurement counts as at
#'   the ceiling (default 0.99).
#' @param bw_check Logical; run the bandwidth check (default \code{TRUE}).
#' @param estimator Kernel estimator settings: a list with elements \code{bw},
#'   \code{engine}, \code{eval_n}, and \code{loo} as returned by
#'   \code{recommended_estimator()}. Defaults to the settings the study used at
#'   this dimensionality, \code{recommended_estimator(length(features))}.
#' @param min_tokens Minimum tokens in the smaller category for a speaker to be
#'   measured at all (default 10). Speakers below it, or without exactly two
#'   observed categories, are left out with a message; the licensing floors
#'   above are applied to the speakers that remain.
#' @param eval_seed Optional integer seed used when \code{estimator$eval_n}
#'   subsamples the evaluation points; makes the estimates reproducible.
#' @param chunk_size Chunk size for \code{engine = "fast_diag"}.
#'
#' @return A tibble of class \code{"phontrast_ranking"}, one row per measured
#'   speaker, sorted by \code{sqrt_jsd} (largest first), with columns
#'   \code{group}, \code{n_tokens}, \code{n_min}, \code{sqrt_jsd},
#'   \code{pillai}, \code{shared_mass}, \code{at_ceiling}, \code{pr_jsd},
#'   \code{pr_pillai}, \code{rank_diff} (\code{pr_jsd - pr_pillai}),
#'   \code{flag}, \code{rank_licensed}, \code{flag_licensed},
#'   \code{rank_basis}, and, with \code{bw_check = TRUE},
#'   \code{sqrt_jsd_half}, \code{sqrt_jsd_double}, \code{bw_shift},
#'   \code{sign_change}, and \code{set_aside}. The protocol settings, the
#'   floors, and the cleaned tokens the ranking was computed from are attached
#'   as \code{attr(x, "protocol")} (so \code{inspect_contrast()} can redraw a
#'   speaker without the original data); \code{print()}
#'   summarizes them. \code{plot()} and \code{ggplot2::autoplot()} draw the
#'   rank-agreement plot.
#'
#' @seealso \code{percentile_rank()}, \code{protocol_floors()},
#'   \code{recommended_estimator()}, \code{inspect_contrast()},
#'   \code{plot_contrast()}, \code{phontrast()}.
#'
#' @examples
#' set.seed(2026)
#' gaps <- c(30, 60, 90, 120, 150, 180, 210, 240)
#' cohort <- do.call(rbind, lapply(seq_along(gaps), function(i) {
#'   data.frame(
#'     speaker = sprintf("s%02d", i),
#'     vowel = rep(c("ih", "eh"), each = 50),
#'     f1 = c(rnorm(50, 500, 55), rnorm(50, 500 + gaps[i], 55)),
#'     f2 = c(rnorm(50, 1900, 120), rnorm(50, 1900 - gaps[i], 120))
#'   )
#' }))
#' ranking <- rank_contrasts(cohort, c("f1", "f2"), "vowel", "speaker",
#'                           bw_check = FALSE)
#' ranking
#' @export
rank_contrasts <- function(data,
                           features,
                           category_col,
                           group_col,
                           margin = 0.25,
                           ceiling = 0.99,
                           bw_check = TRUE,
                           estimator = recommended_estimator(length(features)),
                           min_tokens = 10,
                           eval_seed = NULL,
                           chunk_size = 1000L) {
  if (missing(group_col) || is.null(group_col)) {
    stop(
      "`group_col` is required: rank_contrasts() ranks the speakers (or other ",
      "units) named by `group_col`.",
      call. = FALSE
    )
  }
  .validate_metric_inputs(data, features, category_col, group_col)
  group_col <- .check_group_cols(group_col)
  .check_bool(bw_check, "bw_check")
  .check_positive_count(min_tokens, "min_tokens")
  .check_unit_fraction(margin, "margin")
  .check_unit_fraction(ceiling, "ceiling")
  .check_positive_count(chunk_size, "chunk_size")
  est <- .check_estimator_spec(estimator)
  d <- length(features)
  floors <- protocol_floors(d)
  min_per_category <- max(min_tokens, .kde_min_category_tokens(d))

  df <- .metric_data(data, c(group_col, category_col, features))
  .check_numeric_features(df, features)
  groups <- .split_groups(df, group_col)

  rows <- lapply(groups, function(df_g) {
    .rank_contrasts_one(
      df_g = df_g, features = features, category_col = category_col,
      group_col = group_col, est = est, bw_check = bw_check,
      min_per_category = min_per_category, eval_seed = eval_seed,
      chunk_size = chunk_size
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  dropped <- out[!out$measured, , drop = FALSE]
  out <- out[out$measured, , drop = FALSE]
  if (nrow(dropped)) {
    message(
      "rank_contrasts(): ", nrow(dropped), " group(s) with fewer than ",
      min_per_category, " tokens in a category, or without exactly two observed ",
      "categories, were not measured: ", paste(dropped$group, collapse = ", "), "."
    )
  }
  if (nrow(out) < 2L) {
    stop(
      "rank_contrasts() needs at least two measured groups to rank; got ",
      nrow(out), ". One speaker measured once on one contrast gives nothing to rank.",
      call. = FALSE
    )
  }
  n_failed <- sum(is.na(out$sqrt_jsd) | is.na(out$pillai))
  if (n_failed) {
    warning(
      "rank_contrasts(): ", n_failed, " of ", nrow(out),
      " group(s) could not be fully estimated and carry NA.",
      call. = FALSE
    )
  }

  # ---- Step 2: percentile ranks over the speakers below the ceiling --------
  out$at_ceiling <- !is.na(out$sqrt_jsd) & out$sqrt_jsd >= ceiling
  out$pr_jsd <- percentile_rank(out$sqrt_jsd, exclude = out$at_ceiling)
  out$pr_pillai <- percentile_rank(out$pillai, exclude = out$at_ceiling)
  out$rank_diff <- out$pr_jsd - out$pr_pillai
  k <- sum(!is.na(out$pr_jsd))

  # ---- Licensing ------------------------------------------------------------
  out$rank_licensed <- out$n_min >= floors$rank_floor
  out$flag_licensed <- out$n_min >= floors$flag_floor & .flag_defined(d, out$n_min)
  out$rank_basis <- ifelse(out$rank_licensed, "sqrt_jsd", "pillai")

  # ---- Step 3: the agreement flag, read only where licensed ----------------
  out$flag <- ifelse(out$flag_licensed, abs(out$rank_diff) >= margin, NA)

  # ---- Bandwidth check --------------------------------------------------------
  if (bw_check) {
    pr_half <- percentile_rank(out$sqrt_jsd_half, exclude = out$at_ceiling)
    pr_double <- percentile_rank(out$sqrt_jsd_double, exclude = out$at_ceiling)
    out$bw_shift <- pr_half - pr_double
    diff_half <- pr_half - out$pr_pillai
    diff_double <- pr_double - out$pr_pillai
    out$sign_change <- out$flag %in% TRUE & (diff_half * diff_double < 0) %in% TRUE
    out$set_aside <- abs(out$bw_shift) >= margin | out$sign_change
  } else {
    out$sqrt_jsd_half <- NULL
    out$sqrt_jsd_double <- NULL
  }

  if (k < 8L) {
    warning(
      "rank_contrasts(): only ", k, " speakers enter the ranking; eight is the ",
      "smallest useful set (the calibration used sets of 32).",
      call. = FALSE
    )
  }

  out$measured <- NULL
  out <- out[order(-out$sqrt_jsd, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  first <- c("group", "n_tokens", "n_min", "sqrt_jsd", "pillai", "shared_mass",
             "at_ceiling", "pr_jsd", "pr_pillai", "rank_diff", "flag",
             "rank_licensed", "flag_licensed", "rank_basis")
  out <- out[, c(first, setdiff(names(out), first)), drop = FALSE]
  out <- tibble::as_tibble(out)

  attr(out, "protocol") <- list(
    features = features, d = d, category_col = category_col,
    group_col = group_col, k = k, n_ceiling = sum(out$at_ceiling),
    n_dropped = nrow(dropped), margin = margin, ceiling = ceiling,
    estimator = est, bw_check = bw_check, floors = floors,
    eval_seed = eval_seed, data = df
  )
  class(out) <- unique(c("phontrast_ranking", class(out)))
  out
}

# One speaker: sqrt(JSD) and shared mass from a single kernel pass, Pillai on
# the same tokens, and the two bandwidth-check re-estimates.
.rank_contrasts_one <- function(df_g, features, category_col, group_col, est,
                                bw_check, min_per_category, eval_seed,
                                chunk_size) {
  counts <- .observed_category_counts(df_g[[category_col]])
  n_min <- if (length(counts)) min(as.integer(counts)) else 0L
  row <- data.frame(
    group = .group_label(df_g, group_col),
    n_tokens = nrow(df_g),
    n_min = n_min,
    sqrt_jsd = NA_real_,
    pillai = NA_real_,
    shared_mass = NA_real_,
    sqrt_jsd_half = NA_real_,
    sqrt_jsd_double = NA_real_,
    measured = FALSE,
    stringsAsFactors = FALSE
  )
  if (length(counts) != 2L || n_min < min_per_category) {
    return(row)
  }
  row$measured <- TRUE

  kernel_sqrt_jsd <- function(bw, bw_scale, with_overlap = FALSE) {
    mc <- tryCatch(
      .kde_mc_pair(
        data = df_g, features = features, category_col = category_col,
        bw = bw, eval_n = est$eval_n, eval_seed = eval_seed,
        engine = est$engine, chunk_size = chunk_size,
        metric = "rank_contrasts()", bw_scale = bw_scale
      ),
      error = function(e) NULL
    )
    if (is.null(mc)) {
      return(list(sqrt_jsd = NA_real_, shared_mass = NA_real_))
    }
    list(
      sqrt_jsd = sqrt(.jsd_mc(mc, loo = est$loo)),
      shared_mass = if (with_overlap) .overlap_mc(mc) else NA_real_
    )
  }

  main <- kernel_sqrt_jsd(est$bw, 1, with_overlap = TRUE)
  row$sqrt_jsd <- main$sqrt_jsd
  row$shared_mass <- main$shared_mass
  row$pillai <- tryCatch(
    pillai_overlap(df_g, features, category_col)$pillai,
    error = function(e) NA_real_
  )
  if (bw_check) {
    row$sqrt_jsd_half <- kernel_sqrt_jsd("scott.diag", 0.5)$sqrt_jsd
    row$sqrt_jsd_double <- kernel_sqrt_jsd("scott.diag", 2)$sqrt_jsd
  }
  row
}

# The agreement test is undefined when dimensionality outruns the sample
# (d >= 2n): Pillai's within-class error matrix has no full rank left.
.flag_defined <- function(d, n_min) {
  d < 2 * n_min
}

.check_unit_fraction <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x > 1) {
    stop("`", arg, "` must be a single number in (0, 1].", call. = FALSE)
  }
  invisible(x)
}

.check_estimator_spec <- function(estimator) {
  needed <- c("bw", "engine", "eval_n", "loo")
  if (!is.list(estimator) || !all(needed %in% names(estimator))) {
    stop(
      "`estimator` must be a list with elements `bw`, `engine`, `eval_n`, and ",
      "`loo`, as returned by recommended_estimator().",
      call. = FALSE
    )
  }
  bw <- estimator$bw
  if (!is.character(bw) || length(bw) != 1L ||
      !bw %in% c("Hpi", "Hscv", "Hpi.diag", "scott.diag")) {
    stop(
      "`estimator$bw` must be one of \"Hpi\", \"Hscv\", \"Hpi.diag\", \"scott.diag\".",
      call. = FALSE
    )
  }
  engine <- .match_kde_engine(estimator$engine)
  if (!is.null(estimator$eval_n)) {
    .check_positive_count(estimator$eval_n, "estimator$eval_n")
  }
  .check_bool(estimator$loo, "estimator$loo")
  list(bw = bw, engine = engine, eval_n = estimator$eval_n, loo = estimator$loo)
}

#' Sample-size floors of the ranking protocol
#'
#' The token floors per category at which the simulation study behind
#' \code{rank_contrasts()} licenses a \eqn{\sqrt{JSD}} ranking and the
#' Pillai--\eqn{\sqrt{JSD}} agreement flag, by number of features \code{d}
#' (Berry, under review, Secs. V.D and VII.A). Dimensionalities the study did
#' not calibrate take the floors of the next higher calibrated dimensionality.
#'
#' @param d Positive integer; the number of features.
#'
#' @return A list with \code{d}; \code{rank_floor}, the tokens per category
#'   from which ranking by \eqn{\sqrt{JSD}} is licensed (50 at \code{d <= 2},
#'   200 at \code{d = 3, 4}, 500 at \code{d = 5} to \code{8}, \code{Inf}
#'   above eight dimensions, where the study gives ordering evidence only);
#'   \code{flag_floor}, the tokens per category from which the agreement flag
#'   may be read (100 at \code{d <= 2}, 200 at \code{d = 3, 4}, \code{Inf}
#'   from \code{d = 5}, since from eight dimensions the test returns agreement
#'   whatever the data contain); \code{calibrated_at}, the calibrated
#'   dimensionality the floors come from; and \code{ordering_only}, \code{TRUE}
#'   above eight dimensions.
#'
#' @examples
#' protocol_floors(2)
#' protocol_floors(4)$flag_floor
#' @export
protocol_floors <- function(d) {
  if (!is.numeric(d) || length(d) != 1L || !is.finite(d) || d < 1 || d != round(d)) {
    stop("`d` must be a single positive integer (the number of features).", call. = FALSE)
  }
  d <- as.integer(d)
  if (d <= 2L) {
    list(d = d, rank_floor = 50, flag_floor = 100, calibrated_at = 2L,
         ordering_only = FALSE)
  } else if (d <= 4L) {
    list(d = d, rank_floor = 200, flag_floor = 200, calibrated_at = 4L,
         ordering_only = FALSE)
  } else if (d <= 8L) {
    list(d = d, rank_floor = 500, flag_floor = Inf, calibrated_at = 8L,
         ordering_only = FALSE)
  } else {
    list(d = d, rank_floor = Inf, flag_floor = Inf, calibrated_at = NA_integer_,
         ordering_only = TRUE)
  }
}

#' @export
print.phontrast_ranking <- function(x, ...) {
  p <- attr(x, "protocol")
  if (is.null(p)) {
    return(NextMethod())
  }
  fmt_floor <- function(f) if (is.finite(f)) format(f) else "none"
  est <- p$estimator
  cat(sprintf(
    "<phontrast ranking: %d speaker%s, %d feature%s (%s)>\n",
    nrow(x), if (nrow(x) == 1L) "" else "s", p$d, if (p$d == 1L) "" else "s",
    paste(p$features, collapse = ", ")
  ))
  cat(sprintf(
    "sqrt(JSD), Pillai, and shared mass from the same tokens; kernel: %s / %s, %s, leave-one-out %s.\n",
    est$bw, est$engine,
    if (is.null(est$eval_n)) "all tokens" else paste(est$eval_n, "evaluation tokens"),
    if (isTRUE(est$loo)) "on" else "off"
  ))
  cat(sprintf(
    "Ranked: %d%s. Floors at d = %d: rank from %s, flag from %s tokens per category%s.\n",
    p$k,
    if (p$n_ceiling) sprintf(" (%d at ceiling, sqrt(JSD) >= %s)", p$n_ceiling, format(p$ceiling)) else "",
    p$d, fmt_floor(p$floors$rank_floor), fmt_floor(p$floors$flag_floor),
    if (isTRUE(p$floors$ordering_only)) "; above eight dimensions: ordering evidence only" else ""
  ))
  below_rank <- sum(!x$rank_licensed)
  unread <- sum(is.na(x$flag) & !x$at_ceiling)
  flagged <- x$group[x$flag %in% TRUE]
  cat(sprintf(
    "Below rank floor (rank with Pillai): %d. Flag not readable: %d. Flagged (|rank_diff| >= %s): %d%s.\n",
    below_rank, unread, format(p$margin), length(flagged),
    if (length(flagged)) paste0(" -> ", paste(flagged, collapse = ", ")) else ""
  ))
  if (isTRUE(p$bw_check)) {
    aside <- x$group[x$set_aside %in% TRUE]
    cat(sprintf(
      "Bandwidth check (scott.diag x0.5 / x2): %d set aside%s.\n",
      length(aside),
      if (length(aside)) paste0(" -> ", paste(aside, collapse = ", ")) else ""
    ))
  } else {
    cat("Bandwidth check not run (bw_check = FALSE).\n")
  }
  if (p$n_dropped) {
    cat(sprintf("Not measured (too few tokens or not two categories): %d.\n", p$n_dropped))
  }
  cat("\n")
  y <- x
  class(y) <- setdiff(class(y), "phontrast_ranking")
  attr(y, "protocol") <- NULL
  print(y, ...)
  invisible(x)
}
