# Shared test fixtures for viewseurat testsuite.
#
# This file is sourced automatically by testthat before each test file.
# It provides:
#   - make_counts_matrix()    simple dense matrix with dimnames
#   - make_mock_v4_assay()    minimal S4 object simulating a v4-style Assay
#                             (no 'layers' slot, direct 'counts'/'data' slots)

#' Build a minimal counts matrix suitable for CreateSeuratObject()
#' Returns a sparse dgCMatrix to avoid SeuratObject v5 coercion warnings.
make_counts_matrix <- function(nrow = 2L, ncol = 3L) {
  mat <- matrix(
    seq_len(nrow * ncol),
    nrow  = nrow,
    dimnames = list(
      paste0("gene", seq_len(nrow)),
      paste0("cell", seq_len(ncol))
    )
  )
  Matrix::Matrix(mat, sparse = TRUE)
}

# Minimal S4 class that mimics a v4-style Assay:
#   - has 'counts' and 'data' as direct slots (no 'layers' slot)
#   - is NOT an Assay5 instance
# This exercises the slot-fallback branch in get_layer_dim() and
# get_assay_data_safe().
if (!isClass("MockV4LikeAssay")) {
  methods::setClass(
    "MockV4LikeAssay",
    slots = c(counts = "matrix", data = "matrix")
  )
}

#' Create a MockV4LikeAssay instance
#' Uses a plain dense matrix (not the sparse fixture) because MockV4LikeAssay
#' slots are typed as "matrix" to mirror actual v4 Assay slot types.
make_mock_v4_assay <- function(nrow = 2L, ncol = 3L) {
  mat <- matrix(
    seq_len(nrow * ncol),
    nrow = nrow,
    dimnames = list(
      paste0("gene", seq_len(nrow)),
      paste0("cell", seq_len(ncol))
    )
  )
  methods::new("MockV4LikeAssay", counts = mat, data = mat)
}
