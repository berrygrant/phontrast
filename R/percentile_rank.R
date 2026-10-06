#' Percentile ranks with averaged ties
#'
#' Converts a vector of measurements into percentile ranks,
#' \eqn{(\bar{r} - 1/2) / k}, where \eqn{\bar{r}} is the average rank (ties
#' share the mean of the ranks they span) and \eqn{k} is the number of
#' measurements ranked. This is the rank scale on which the agreement check of
#' \code{rank_contrasts()} compares Jensen-Shannon distance with Pillai:
#' a difference of 0.25 means the two measures place a speaker a quarter of the
#' ordering apart. Ranks are ascending, so larger values of \code{x} receive
#' larger percentile ranks. Measurements that are \code{NA} or flagged in
#' \code{exclude} are left out of the ranking (they do not count toward
#' \eqn{k}) and come back as \code{NA}.
#'
#' @param x Numeric vector of measurements to rank.
#' @param exclude Optional logical vector the same length as \code{x}:
#'   \code{TRUE} removes a measurement from the ranking before the ranks are
#'   computed. For example, \code{exclude = x >= 0.99} drops Jensen-Shannon
#'   distances at the ceiling. \code{NA} entries are treated as \code{FALSE}.
#'
#' @return Numeric vector the same length as \code{x} with percentile ranks in
#'   \eqn{(0, 1)} for the ranked measurements and \code{NA} elsewhere.
#'
#' @examples
#' percentile_rank(c(0.2, 0.5, 0.5, 0.9))
#' percentile_rank(c(0.2, 0.5, 0.5, 0.995), exclude = c(0.2, 0.5, 0.5, 0.995) >= 0.99)
#' @export
percentile_rank <- function(x, exclude = NULL) {
  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (is.null(exclude)) {
    exclude <- rep(FALSE, length(x))
  }
  if (!is.logical(exclude) || length(exclude) != length(x)) {
    stop("`exclude` must be a logical vector the same length as `x`.", call. = FALSE)
  }
  exclude[is.na(exclude)] <- FALSE

  keep <- !exclude & !is.na(x)
  out <- rep(NA_real_, length(x))
  k <- sum(keep)
  if (k == 0L) {
    return(out)
  }
  r <- rank(x[keep], ties.method = "average")
  out[keep] <- (r - 0.5) / k
  out
}
