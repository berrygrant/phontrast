test_that("citation metadata matches the installed package release", {
  package_version <- as.character(utils::packageVersion("phontrast"))
  citation_entry <- utils::citation("phontrast")[[1L]]

  expect_identical(unname(citation_entry$version), package_version)
  # The text rendering wraps lines, and a long (development) version string can
  # wrap between "version" and the number, so compare with whitespace collapsed.
  expect_match(
    gsub("\\s+", " ", format(citation_entry, style = "text")),
    paste("version", package_version),
    fixed = TRUE
  )

  expected_dois <- c(
    "2.3.1" = "10.5281/zenodo.21795954",
    "2.4.0" = "10.5281/zenodo.21864533",
    "2.5.0" = "10.5281/zenodo.23193095"
  )
  expected_doi <- unname(expected_dois[package_version])
  if (is.na(expected_doi)) {
    expected_doi <- "10.5281/zenodo.20816585"
  }
  expect_identical(unname(citation_entry$doi), expected_doi)
})
