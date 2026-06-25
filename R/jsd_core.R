#' Kullback–Leibler divergence for discrete distributions
#'
#' Computes KL(p || q) in bits for discrete probability vectors.
#' Zero-probability events in `p` contribute zero; positive mass in `p`
#' where `q` is zero returns `Inf`.
#'
#' @param p,q Numeric probability vectors of the same length.
#'
#' @return A single numeric value: the KL divergence in bits.
#' @export
kl_div <- function(p, q) {
  if (length(p) != length(q)) {
    stop("`p` and `q` must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(p)) || any(!is.finite(q)) || any(p < 0) || any(q < 0)) {
    stop("`p` and `q` must be finite, non-negative vectors.", call. = FALSE)
  }

  keep <- p > 0
  if (any(q[keep] == 0)) {
    return(Inf)
  }
  sum(p[keep] * log2(p[keep] / q[keep]))
}

#' Jensen–Shannon divergence for discrete distributions
#'
#' Computes JSD(p, q) in bits for discrete probability vectors.
#' The value is bounded in \code{[0, 1]} for equally weighted mixtures.
#'
#' @inheritParams kl_div
#'
#' @return A single numeric value: the Jensen–Shannon divergence in bits.
#' @export
jsd <- function(p, q) {
  if (length(p) != length(q)) {
    stop("`p` and `q` must have the same length.", call. = FALSE)
  }
  if (any(!is.finite(p)) || any(!is.finite(q)) || any(p < 0) || any(q < 0)) {
    stop("`p` and `q` must be finite, non-negative vectors.", call. = FALSE)
  }
  if (sum(p) <= 0 || sum(q) <= 0) {
    stop("`p` and `q` must each have positive mass.", call. = FALSE)
  }

  p <- p / sum(p)
  q <- q / sum(q)
  m <- 0.5 * (p + q)
  0.5 * kl_div(p, m) + 0.5 * kl_div(q, m)
}
