#' Extract MFCCs for vowel segments
#'
#' Adds MFCC feature columns (e.g., mfcc1..mfcc13) to a data frame that
#' contains audio file paths and (optionally) segment boundaries.
#'
#' This function uses \pkg{tuneR} to read WAV files and \pkg{seewave} to
#' compute MFCCs. Both packages are optional; if not installed, an
#' informative error is raised.
#'
#' @param data Data frame containing audio paths and segment boundaries.
#' @param file_col String; column name containing WAV file paths.
#' @param start_col Optional string; column name containing segment start time
#'   (in seconds). If \code{NULL}, the full file is used.
#' @param end_col Optional string; column name containing segment end time
#'   (in seconds). If \code{NULL}, the full file is used.
#' @param fs Optional numeric; override sampling rate (Hz). If \code{NULL},
#'   uses the file's sampling rate.
#' @param numcep Integer; number of MFCC coefficients to return.
#' @param prefix String; prefix for output columns (default \code{"mfcc"}).
#' @param ... Additional arguments passed to \code{seewave::mfcc()}.
#'
#' @return The input data frame with added MFCC columns.
#' @export
extract_mfcc <- function(data,
                         file_col,
                         start_col = NULL,
                         end_col   = NULL,
                         fs        = NULL,
                         numcep    = 13,
                         prefix    = "mfcc",
                         ...) {

  if (!requireNamespace("tuneR", quietly = TRUE)) {
    stop("extract_mfcc(): package 'tuneR' is required but not installed.")
  }
  if (!requireNamespace("seewave", quietly = TRUE)) {
    stop("extract_mfcc(): package 'seewave' is required but not installed.")
  }
  if (!file_col %in% names(data)) {
    stop("extract_mfcc(): `file_col` must be a column in `data`.")
  }
  if (!is.null(start_col) && !start_col %in% names(data)) {
    stop("extract_mfcc(): `start_col` must be a column in `data`.")
  }
  if (!is.null(end_col) && !end_col %in% names(data)) {
    stop("extract_mfcc(): `end_col` must be a column in `data`.")
  }

  n <- nrow(data)
  out <- matrix(NA_real_, nrow = n, ncol = numcep)
  colnames(out) <- paste0(prefix, seq_len(numcep))

  for (i in seq_len(n)) {
    path <- data[[file_col]][i]
    if (is.na(path) || !nzchar(path)) {
      next
    }

    wave <- tryCatch(
      tuneR::readWave(path),
      error = function(e) NULL
    )
    if (is.null(wave)) {
      next
    }

    # Extract segment if boundaries are provided
    if (!is.null(start_col) && !is.null(end_col)) {
      start_t <- data[[start_col]][i]
      end_t   <- data[[end_col]][i]
      if (is.finite(start_t) && is.finite(end_t) && end_t > start_t) {
        wave <- tryCatch(
          tuneR::extractWave(wave, from = start_t, to = end_t, xunit = "time"),
          error = function(e) wave
        )
      }
    }

    # Convert to mono if needed
    if (isTRUE(wave@stereo)) {
      wave <- tuneR::mono(wave, which = "left")
    }

    fs_use <- if (is.null(fs)) wave@samp.rate else fs
    x <- wave@left

    mf <- tryCatch(
      seewave::mfcc(x = x, f = fs_use, ...),
      error = function(e) NULL
    )
    if (is.null(mf)) {
      next
    }

    mf <- as.matrix(mf)
    vec <- rep(NA_real_, numcep)

    # Try to interpret dimensions (frames x coeffs is typical)
    if (ncol(mf) >= numcep) {
      mf_use <- mf[, seq_len(numcep), drop = FALSE]
      vec <- colMeans(mf_use)
    } else if (nrow(mf) >= numcep) {
      mf_use <- mf[seq_len(numcep), , drop = FALSE]
      vec <- rowMeans(mf_use)
    } else if (ncol(mf) > 1) {
      v <- colMeans(mf)
      vec[seq_len(length(v))] <- v
    } else if (nrow(mf) > 1) {
      v <- rowMeans(mf)
      vec[seq_len(length(v))] <- v
    }

    out[i, ] <- vec
  }

  cbind(data, out)
}
