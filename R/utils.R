.check_columns <- function(data, cols, arg = "columns") {
  missing <- setdiff(cols, names(data))
  if (length(missing)) {
    stop(
      "`", arg, "` must exist in `data`: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

.metric_data <- function(data, cols) {
  .check_columns(data, cols)
  data[stats::complete.cases(data[, cols, drop = FALSE]), cols, drop = FALSE]
}

.two_levels <- function(x, arg = "category") {
  levs <- unique(x)
  levs <- levs[!is.na(levs)]
  if (length(levs) != 2L) {
    stop("`", arg, "` must have exactly two non-missing values.", call. = FALSE)
  }
  levs
}

.empty_group_jsd <- function() {
  tibble::tibble(group = character(), n_tokens = integer(), jsd = numeric())
}

.empty_group_bhatt <- function() {
  data.frame(
    group = character(),
    n_tokens = integer(),
    bhatt_dist = numeric(),
    bhatt_affinity = numeric(),
    stringsAsFactors = FALSE
  )
}
