# Tests for get_layer_dim() --------------------------------------------------
# The function handles two assay layouts:
#   v5 / Assay5 – matrix accessed via the @layers list
#   v4-style    – matrix stored in a direct named slot (counts, data, …)
# SCTAssay is v4-style and is covered by the mock path below.

# -- v5 Assay5 path -----------------------------------------------------------

test_that("get_layer_dim returns correct dims for Assay5 counts layer", {
  skip_if_not_installed("SeuratObject")

  assay5 <- SeuratObject::CreateAssay5Object(
    counts = make_counts_matrix(nrow = 3L, ncol = 4L)
  )
  expect_equal(get_layer_dim(assay5, "counts"), c(3L, 4L))
})

test_that("get_layer_dim returns c(0L, 0L) for absent layer in Assay5", {
  skip_if_not_installed("SeuratObject")

  assay5 <- SeuratObject::CreateAssay5Object(
    counts = make_counts_matrix()
  )
  expect_equal(get_layer_dim(assay5, "nonexistent"), c(0L, 0L))
})

test_that("get_layer_dim returns c(0L, 0L) for unpopulated scale.data in Assay5", {
  skip_if_not_installed("SeuratObject")

  # scale.data is not populated by CreateAssay5Object
  assay5 <- SeuratObject::CreateAssay5Object(
    counts = make_counts_matrix()
  )
  expect_equal(get_layer_dim(assay5, "scale.data"), c(0L, 0L))
})

# -- v4-style / slot-based path (MockV4LikeAssay from helper-fixtures.R) ------

test_that("get_layer_dim returns correct dims for v4-style slot-based assay", {
  mock <- make_mock_v4_assay(nrow = 2L, ncol = 3L)
  expect_equal(get_layer_dim(mock, "counts"), c(2L, 3L))
  expect_equal(get_layer_dim(mock, "data"),   c(2L, 3L))
})

test_that("get_layer_dim returns c(0L, 0L) for slot absent in v4-style assay", {
  mock <- make_mock_v4_assay()
  # MockV4LikeAssay has no 'scale.data' slot
  expect_equal(get_layer_dim(mock, "scale.data"), c(0L, 0L))
  expect_equal(get_layer_dim(mock, "nonexistent"), c(0L, 0L))
})
