test_that("validate_seurat_object rejects non-Seurat objects", {
  expect_error(validate_seurat_object(list()),       "not a Seurat object")
  expect_error(validate_seurat_object(NULL),         "not a Seurat object")
  expect_error(validate_seurat_object(42L),          "not a Seurat object")
  expect_error(validate_seurat_object("string"),     "not a Seurat object")
  expect_error(validate_seurat_object(data.frame()), "not a Seurat object")
})

test_that("validate_seurat_object returns TRUE for a valid Seurat object", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())
  expect_true(validate_seurat_object(obj))
})

test_that("validate_seurat_object errors when the object has no cells", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())

  # SeuratObject v5 rejects subsetting to 0 cells at the [ operator level,
  # so we can't construct obj_empty via obj[, integer(0)].  Skip gracefully
  # if this version enforces that constraint.
  obj_empty <- tryCatch(
    obj[, integer(0)],
    error = function(e) skip(paste("SeuratObject prevents empty subset:", conditionMessage(e)))
  )

  expect_error(validate_seurat_object(obj_empty), "no cells")
})

test_that("validate_seurat_object errors when the object has no assays", {
  skip_if_not_installed("SeuratObject")

  obj <- SeuratObject::CreateSeuratObject(counts = make_counts_matrix())
  obj@assays <- list()

  expect_error(validate_seurat_object(obj), "no assays")
})
