# Tests that app.R sources all R/ files needed at runtime -----------------------
#
# Regression: metadata_utils.R was missing from app.R's source() list, so
# filter_metadata_columns() would be undefined when running via RStudio
# "Run App" (which uses app.R directly rather than devtools::load_all()).
#
# Files legitimately exempt from app.R:
# - viewseurat-package.R: pure roxygen @import declarations, no runtime code

# ViewSeurat.R: public entry point called by users; app.R creates the app
#   directly and never calls ViewSeurat()
# viewseurat-package.R: pure roxygen @import declarations, no runtime code
app_source_exemptions <- c("ViewSeurat.R", "viewseurat-package.R")

test_that("app.R sources every R/ file that defines runtime functions", {
  app_lines <- readLines("../../app.R")

  sourced_files <- regmatches(
    app_lines,
    regexpr('(?<=source\\(")R/[^"]+', app_lines, perl = TRUE)
  )
  sourced_bases <- basename(sourced_files)

  r_files <- list.files("../../R", pattern = "\\.R$")
  expected <- setdiff(r_files, app_source_exemptions)

  missing <- setdiff(expected, sourced_bases)
  expect_equal(
    missing, character(0),
    info = paste("R/ files not sourced in app.R:", paste(missing, collapse = ", "))
  )
})
