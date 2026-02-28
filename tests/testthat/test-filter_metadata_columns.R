test_that("filter_metadata_columns returns all columns when search is empty", {
  df <- data.frame(nCount_RNA = 1:3, nFeature_RNA = 4:6, percent.mt = 7:9,
                   row.names = c("AACG", "TTCG", "GGCA"))
  expect_equal(filter_metadata_columns(df, ""), df)
  expect_equal(filter_metadata_columns(df, NULL), df)
  expect_equal(filter_metadata_columns(df, "   "), df)
})

test_that("filter_metadata_columns subsets by partial match", {
  df <- data.frame(nCount_RNA = 1:3, nFeature_RNA = 4:6, percent.mt = 7:9,
                   row.names = c("AACG", "TTCG", "GGCA"))
  result <- filter_metadata_columns(df, "RNA")
  expect_equal(ncol(result), 2)
  expect_true(all(c("nCount_RNA", "nFeature_RNA") %in% colnames(result)))
  expect_equal(rownames(result), c("AACG", "TTCG", "GGCA"))
})

test_that("filter_metadata_columns is case-insensitive", {
  df <- data.frame(nCount_RNA = 1:3, percent.mt = 4:6,
                   row.names = c("A", "B", "C"))
  result <- filter_metadata_columns(df, "rna")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "nCount_RNA")
})

test_that("filter_metadata_columns returns zero columns when no match", {
  df <- data.frame(nCount_RNA = 1:3, percent.mt = 4:6,
                   row.names = c("A", "B", "C"))
  result <- filter_metadata_columns(df, "zzz_nonexistent")
  expect_equal(ncol(result), 0)
  expect_equal(nrow(result), 3)
  expect_equal(rownames(result), c("A", "B", "C"))
})

test_that("filter_metadata_columns returns data.frame for single column match", {
  df <- data.frame(nCount_RNA = 1:3, nFeature_RNA = 4:6, percent.mt = 7:9,
                   row.names = c("A", "B", "C"))
  result <- filter_metadata_columns(df, "percent")
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "percent.mt")
})

test_that("filter_metadata_columns works with single-row data.frame", {
  df <- data.frame(nCount_RNA = 10, percent.mt = 0.5,
                   row.names = "AACG")
  result <- filter_metadata_columns(df, "RNA")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 1)
  expect_equal(rownames(result), "AACG")
})

test_that("filter_metadata_columns uses literal matching not regex", {
  df <- data.frame(percent.mt = 1:3, percent_mt = 4:6, percentXmt = 7:9,
                   row.names = c("A", "B", "C"))
  result <- filter_metadata_columns(df, "percent.mt")
  expect_equal(ncol(result), 1)
  expect_equal(colnames(result), "percent.mt")
})
