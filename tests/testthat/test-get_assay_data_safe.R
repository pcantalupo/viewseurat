# Tests for get_assay_data_safe() --------------------------------------------
# The function wraps LayerData for Assay5 objects and falls back to direct
# slot access for v4-style assays (no @layers slot, not Assay5).

# -- Assay5 (v5) path ---------------------------------------------------------

test_that("get_assay_data_safe returns a matrix for a valid assay/layer", {
  skip_if_not_installed("SeuratObject")

  obj    <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix(3L, 4L))
  result <- get_assay_data_safe(obj, "RNA", "counts")

  expect_false(is.null(result))
  expect_true(is.matrix(result) || inherits(result, "Matrix"))
  expect_equal(dim(result), c(3L, 4L))
})

test_that("get_assay_data_safe returns NULL for a nonexistent assay name", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())
  expect_null(get_assay_data_safe(obj, "MISSING_ASSAY", "counts"))
})

test_that("get_assay_data_safe returns NULL for a nonexistent layer name", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())
  expect_null(get_assay_data_safe(obj, "RNA", "nonexistent_layer"))
})

test_that("get_assay_data_safe subsets rows with max_features", {
  skip_if_not_installed("SeuratObject")

  obj    <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix(5L, 4L))
  result <- get_assay_data_safe(obj, "RNA", "counts", max_features = 2L)

  expect_equal(nrow(result), 2L)
  expect_equal(ncol(result), 4L)
})

test_that("get_assay_data_safe subsets columns with max_cells", {
  skip_if_not_installed("SeuratObject")

  obj    <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix(5L, 6L))
  result <- get_assay_data_safe(obj, "RNA", "counts", max_cells = 3L)

  expect_equal(nrow(result), 5L)
  expect_equal(ncol(result), 3L)
})

test_that("get_assay_data_safe does not subset when dims are within limits", {
  skip_if_not_installed("SeuratObject")

  obj    <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix(3L, 4L))
  result <- get_assay_data_safe(obj, "RNA", "counts",
                                max_features = 100L, max_cells = 100L)
  expect_equal(dim(result), c(3L, 4L))
})

# -- v4-style slot-fallback path ----------------------------------------------
# We attach a MockV4LikeAssay (defined in helper-fixtures.R) directly to the
# Seurat object's @assays list to exercise the slot-fallback branch without
# needing a Seurat v3-era installation.

test_that("get_assay_data_safe slot-fallback returns matrix for present slot", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())

  # Inject a mock v4-style assay directly into the assays list.
  # Direct @-slot replacement is safe here: SeuratObject types @assays as a
  # plain list, so arbitrary values can be stored.  If a future SeuratObject
  # version rejects non-Assay values we skip rather than hard-fail.
  tryCatch(
    obj@assays[["mock_v4"]] <- make_mock_v4_assay(nrow = 2L, ncol = 3L),
    error = function(e) skip(paste("@assays assignment rejected:", conditionMessage(e)))
  )

  result <- get_assay_data_safe(obj, "mock_v4", "counts")
  expect_equal(dim(result), c(2L, 3L))
})

test_that("get_assay_data_safe slot-fallback returns NULL for absent slot", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())

  tryCatch(
    obj@assays[["mock_v4"]] <- make_mock_v4_assay(),
    error = function(e) skip(paste("@assays assignment rejected:", conditionMessage(e)))
  )

  # MockV4LikeAssay has no 'scale.data' slot
  expect_null(get_assay_data_safe(obj, "mock_v4", "scale.data"))
})
