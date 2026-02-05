#' Get Matrix Sample
#'
#' Extract a sample of rows and columns from a matrix for display.
#'
#' @param matrix A matrix (dense or sparse)
#' @param max_rows Maximum number of rows to return
#' @param max_cols Maximum number of columns to return
#' @return A dense matrix sample
#' @keywords internal
#' @export
get_matrix_sample <- function(matrix, max_rows = 10, max_cols = 20) {
  # Handle empty matrices - 1:0 in R returns c(1, 0), not empty vector
 nr <- nrow(matrix)
  nc <- ncol(matrix)

  if (nr == 0 || nc == 0) {
    if (inherits(matrix, "dgCMatrix") || inherits(matrix, "sparseMatrix")) {
      return(as.matrix(matrix))
    } else {
      return(matrix)
    }
  }

  row_idx <- seq_len(min(nr, max_rows))
  col_idx <- seq_len(min(nc, max_cols))

  if (inherits(matrix, "dgCMatrix") || inherits(matrix, "sparseMatrix")) {
    return(as.matrix(matrix[row_idx, col_idx, drop = FALSE]))
  } else {
    return(matrix[row_idx, col_idx, drop = FALSE])
  }
}

#' Get Sparsity Information
#'
#' Calculate sparsity statistics for a matrix.
#'
#' @param matrix A matrix (dense or sparse)
#' @return A list with sparsity statistics
#' @keywords internal
#' @export
get_sparsity_info <- function(matrix) {
  total_elements <- as.numeric(nrow(matrix)) * as.numeric(ncol(matrix))

  # Handle empty matrices to avoid division by zero
 if (total_elements == 0) {
    return(list(
      total_elements = 0,
      non_zero_elements = 0,
      zero_elements = 0,
      sparsity_percent = 100,
      memory_mb = 0
    ))
  }

  # Only use slot access for dgCMatrix specifically (not all sparseMatrix subclasses)
 if (inherits(matrix, "dgCMatrix")) {
    non_zero <- length(matrix@x)
    sparsity <- 1 - (non_zero / total_elements)

    # Estimate memory from sparse slots directly to avoid expensive traversal
    memory_bytes <- as.numeric(length(matrix@x) * 8 +  # numeric values
                               length(matrix@i) * 4 +  # integer row indices
                               length(matrix@p) * 4)   # integer column pointers
    info <- list(
      total_elements = total_elements,
      non_zero_elements = non_zero,
      zero_elements = total_elements - non_zero,
      sparsity_percent = round(sparsity * 100, 2),
      memory_mb = round(memory_bytes / 1024^2, 2)
    )
  } else if (inherits(matrix, "sparseMatrix")) {
    # For other sparse matrix types, use generic approach
   non_zero <- Matrix::nnzero(matrix)
    sparsity <- 1 - (non_zero / total_elements)
    info <- list(
      total_elements = total_elements,
      non_zero_elements = non_zero,
      zero_elements = total_elements - non_zero,
      sparsity_percent = round(sparsity * 100, 2),
      memory_mb = as.numeric(object.size(matrix)) / 1024^2
    )
  } else {
    non_zero <- sum(matrix != 0)
    sparsity <- 1 - (non_zero / total_elements)

    info <- list(
      total_elements = total_elements,
      non_zero_elements = non_zero,
      zero_elements = total_elements - non_zero,
      sparsity_percent = round(sparsity * 100, 2),
      memory_mb = as.numeric(object.size(matrix)) / 1024^2
    )
  }

  return(info)
}

#' Get Top Variable Features
#'
#' Get the top variable features for an assay.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay
#' @param n Number of features to return
#' @return Character vector of feature names
#' @keywords internal
#' @export
#' @importFrom Seurat VariableFeatures DefaultAssay<-
get_top_variable_features <- function(obj, assay_name, n = 20) {
  DefaultAssay(obj) <- assay_name
  var_features <- VariableFeatures(obj)

  if (length(var_features) == 0) {
    return(NULL)
  }

  head(var_features, n)
}

#' Format Number
#'
#' Format large numbers with K/M/B suffixes.
#'
#' @param x A numeric value
#' @return A formatted string
#' @keywords internal
#' @export
format_number <- function(x) {
  if (x >= 1e9) {
    return(paste0(round(x / 1e9, 2), "B"))
  } else if (x >= 1e6) {
    return(paste0(round(x / 1e6, 2), "M"))
  } else if (x >= 1e3) {
    return(paste0(round(x / 1e3, 2), "K"))
  } else {
    return(as.character(x))
  }
}

#' Validate Seurat Object
#'
#' Check that a Seurat object is valid for viewing.
#'
#' @param obj An object to validate
#' @return TRUE if valid, otherwise throws an error
#' @keywords internal
#' @export
validate_seurat_object <- function(obj) {
  # Check type first and return early - subsequent checks would fail with cryptic errors
 if (!inherits(obj, "Seurat")) {
    stop("Object is not a Seurat object")
  }

  errors <- c()

  if (ncol(obj) == 0) {
    errors <- c(errors, "Object contains no cells")
  }

  if (length(obj@assays) == 0) {
    errors <- c(errors, "Object contains no assays")
  }

  if (length(errors) > 0) {
    stop(paste(errors, collapse = "; "))
  }

  return(TRUE)
}

#' Get Assay Data Safely with Subsetting
#'
#' Wrapper around LayerData with error handling and optional subsetting.
#' Subsetting at retrieval time is much faster for large/on-disk matrices
#' because LayerData() with features/cells params avoids loading the full matrix.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay
#' @param layer Name of the layer (counts, data, scale.data)
#' @param max_features Maximum features to retrieve (NULL for all)
#' @param max_cells Maximum cells to retrieve (NULL for all)
#' @return Matrix data or NULL if not available
#' @keywords internal
#' @export
#' @importFrom SeuratObject LayerData
get_assay_data_safe <- function(obj, assay_name, layer,
                                max_features = NULL, max_cells = NULL) {
  tryCatch({
    assay <- obj@assays[[assay_name]]

    # Determine subsetting indices
    features <- NULL
    cells <- NULL

    if (!is.null(max_features)) {
      n_features <- nrow(assay)
      if (max_features < n_features) {
        features <- 1:max_features
      }
    }

    if (!is.null(max_cells)) {
      n_cells <- ncol(assay)
      if (max_cells < n_cells) {
        cells <- 1:max_cells
      }
    }

    # LayerData subsets at retrieval - much faster for large/on-disk data
    SeuratObject::LayerData(assay, layer = layer,
                            features = features, cells = cells)
  }, error = function(e) {
    NULL
  })
}

#' Find Ident Label
#'
#' Find the metadata column name that matches the active Idents.
#'
#' @param seurat A Seurat object
#' @return The name of the matching metadata column, or "Unknown"
#' @keywords internal
#' @export
#' @importFrom Seurat Idents
FindIdentLabel <- function(seurat) {
  # Find the metadata column name that matches Idents
  ident_values <- as.character(Idents(seurat))

  for (col in colnames(seurat@meta.data)) {
    if (identical(as.character(seurat@meta.data[[col]]), ident_values)) {
      return(col)
    }
  }

  return("Unknown")
}

#' Get Layer Dimensions
#'
#' Returns the dimensions of a specific layer within an assay, handling both
#' v5 Assay objects (layers slot) and v4-style direct slots. Works on SCTAssay
#' and other assay subclasses that may lack a layers slot.
#'
#' @param assay_obj A Seurat Assay object
#' @param layer Name of the layer (e.g. "counts", "data", "scale.data")
#' @return Integer vector c(nrow, ncol), or c(0L, 0L) if layer is absent
#' @keywords internal
#' @export
get_layer_dim <- function(assay_obj, layer) {
  slotnames <- slotNames(assay_obj)
  mat <- NULL
  if ("layers" %in% slotnames) {
    mat <- assay_obj@layers[[layer]]
  } else if (layer %in% slotnames) {
    mat <- slot(assay_obj, layer)
  }

  if (is.null(mat)) return(c(0L, 0L))
  dim(mat)
}

#' Get Seurat Object Information
#'
#' Generate a comprehensive summary of a Seurat object.
#'
#' @param seurat A Seurat object
#' @return A list with formatted information about the object
#' @export
#' @importFrom Seurat DefaultAssay VariableFeatures Idents
SeuratInfo <- function(seurat) {
  output <- list()

  # Seurat version
  output$version <- paste("Seurat version:", as.character(seurat@version))

  # Graphs
  if (length(seurat@graphs) > 0) {
    output$graphs <- paste("Graphs:", paste(names(seurat@graphs), collapse = ", "))
  } else {
    output$graphs <- "Graphs: None"
  }

  # Reductions
  if (length(seurat@reductions) > 0) {
    assayused <- sapply(names(seurat@reductions), function(name) {
      return(seurat[[name]]@assay.used)
    })
    output$reductions <- paste("Reductions:",
                               paste(names(assayused), " (", assayused, ")", sep = "", collapse = ", "))
  } else {
    output$reductions <- "Reductions: None"
  }

  # Images
  if (length(seurat@images) > 0) {
    output$images <- paste("Images:", paste(names(seurat@images), collapse = ", "))
  } else {
    output$images <- "Images: None"
  }

  # Idents
  output$ident_label <- paste("Ident label:", FindIdentLabel(seurat))

  tab <- table(Idents(seurat))
  df <- data.frame(Count = as.integer(tab))
  rownames(df) <- rownames(tab)
  output$idents_table <- t(df)

  # Assays
  assays <- names(seurat@assays)
  default_assay <- DefaultAssay(seurat)

  slotinfo <- list()
  slots <- c("counts", "data", "scale.data")

  for (assay in assays) {
    assay_obj <- seurat@assays[[assay]]

    defaultassay <- ifelse(default_assay == assay, "YES", "")

    dims <- sapply(slots, function(slot) {
      d <- get_layer_dim(assay_obj, slot)
      paste0(d[1], "x", d[2])
    })

    hvgs <- length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] <- c(defaultassay, dims, hvgs)
  }

  output$assays_table <- data.frame(do.call(rbind, slotinfo))
  colnames(output$assays_table) <- c("default", slots, "HVGs")

  return(output)
}
