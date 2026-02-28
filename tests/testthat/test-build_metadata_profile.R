# Tests for build_metadata_profile() -----------------------------------------
# This function requires only a plain data.frame, so no Seurat objects needed.

# -- Output structure ---------------------------------------------------------

test_that("build_metadata_profile returns a data.frame with the correct columns", {
  result <- build_metadata_profile(data.frame(x = 1:3))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("Column", "Type", "Distribution", "Stats", "Complete"))
})

test_that("build_metadata_profile returns one row per metadata column", {
  meta <- data.frame(
    a = 1:5,
    b = factor(letters[1:5]),
    c = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    d = letters[1:5],
    stringsAsFactors = FALSE
  )
  result <- build_metadata_profile(meta)
  expect_equal(nrow(result), 4L)
  expect_equal(result$Column, colnames(meta))
})

# -- Numeric column -----------------------------------------------------------

test_that("build_metadata_profile classifies numeric columns correctly", {
  meta   <- data.frame(score = c(1.0, 2.5, 3.0, 4.5, 5.0))
  result <- build_metadata_profile(meta)

  expect_equal(result$Type, "numeric")
  expect_true(grepl("min:",  result$Stats))
  expect_true(grepl("med:",  result$Stats))
  expect_true(grepl("mean:", result$Stats))
  expect_true(grepl("max:",  result$Stats))
})

test_that("build_metadata_profile classifies integer columns correctly", {
  meta   <- data.frame(count = 1L:5L)
  result <- build_metadata_profile(meta)

  expect_equal(result$Type, "integer")
  expect_true(grepl("min:", result$Stats))
})

# -- Logical column -----------------------------------------------------------

test_that("build_metadata_profile classifies logical columns correctly", {
  meta   <- data.frame(flag = c(TRUE, FALSE, TRUE, FALSE, TRUE))
  result <- build_metadata_profile(meta)

  expect_equal(result$Type, "logical")
  expect_true(grepl("TRUE:",  result$Stats))
  expect_true(grepl("FALSE:", result$Stats))
})

test_that("build_metadata_profile reports NA count for logical with NAs", {
  meta   <- data.frame(flag = c(TRUE, FALSE, NA))
  result <- build_metadata_profile(meta)

  expect_true(grepl("NA:", result$Stats))
})

# -- Factor column ------------------------------------------------------------

test_that("build_metadata_profile classifies factor columns correctly", {
  meta   <- data.frame(celltype = factor(c("T", "B", "T", "NK", "T")))
  result <- build_metadata_profile(meta)

  expect_equal(result$Type, "factor")
  expect_true(grepl("unique", result$Stats))
})

test_that("build_metadata_profile reports top max_levels for factors", {
  lvls <- paste0("type", 1:10)
  meta <- data.frame(celltype = factor(rep(lvls, each = 2)))

  result_5  <- build_metadata_profile(meta, max_levels = 5L)
  result_3  <- build_metadata_profile(meta, max_levels = 3L)

  expect_true(grepl("showing top 5", result_5$Stats))
  expect_true(grepl("showing top 3", result_3$Stats))
})

# -- Character column ---------------------------------------------------------

test_that("build_metadata_profile classifies character columns correctly", {
  meta   <- data.frame(sample = c("s1", "s2", "s1", "s3"),
                       stringsAsFactors = FALSE)
  result <- build_metadata_profile(meta)

  expect_equal(result$Type, "character")
  expect_true(grepl("unique", result$Stats))
})

# -- Completeness -------------------------------------------------------------

test_that("build_metadata_profile reports completeness as a percentage string", {
  meta   <- data.frame(x = c(1, NA, 3, NA, 5))   # 3/5 complete = 60 %
  result <- build_metadata_profile(meta)

  expect_equal(result$Complete, "60%")
})

test_that("build_metadata_profile reports 100% completeness when no NAs", {
  meta   <- data.frame(x = 1:4)
  result <- build_metadata_profile(meta)

  expect_equal(result$Complete, "100%")
})

# -- All-NA column ------------------------------------------------------------

test_that("build_metadata_profile handles an all-NA numeric column gracefully", {
  meta   <- data.frame(x = c(NA_real_, NA_real_, NA_real_))
  result <- build_metadata_profile(meta)

  expect_equal(result$Stats,    "all NA")
  expect_equal(result$Complete, "0%")
})

test_that("build_metadata_profile handles an all-NA logical column gracefully", {
  meta   <- data.frame(flag = c(NA, NA, NA))
  result <- build_metadata_profile(meta)

  expect_equal(result$Stats,    "all NA")
  expect_equal(result$Complete, "0%")
})

test_that("build_metadata_profile handles an all-NA character column gracefully", {
  meta   <- data.frame(label = c(NA_character_, NA_character_),
                       stringsAsFactors = FALSE)
  result <- build_metadata_profile(meta)

  expect_equal(result$Stats,    "all NA")
  expect_equal(result$Complete, "0%")
})

# -- HTML output sanity -------------------------------------------------------

test_that("build_metadata_profile Distribution column contains HTML", {
  meta   <- data.frame(x = 1:5)
  result <- build_metadata_profile(meta)

  # The distribution is rendered as inline HTML
  expect_true(grepl("<div", result$Distribution))
})

test_that("build_metadata_profile escapes HTML special chars in level names", {
  meta <- data.frame(
    grp = c("<script>", "<script>", "normal"),
    stringsAsFactors = FALSE
  )
  result <- build_metadata_profile(meta)

  # Raw '<script>' must not appear in Distribution HTML
  expect_false(grepl("<script>", result$Distribution, fixed = TRUE))
  expect_true(grepl("&lt;script&gt;", result$Distribution, fixed = TRUE))
})

test_that("build_metadata_profile escapes HTML special chars in factor level names", {
  meta <- data.frame(
    grp = factor(c("<b>bold</b>", "<b>bold</b>", "normal")),
    stringsAsFactors = FALSE
  )
  result <- build_metadata_profile(meta)

  expect_false(grepl("<b>", result$Distribution, fixed = TRUE))
  expect_true(grepl("&lt;b&gt;", result$Distribution, fixed = TRUE))
})
