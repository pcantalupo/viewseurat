#' Get Assay Information
#'
#' Extract information about a specific assay from a Seurat object.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay to get information for
#' @return A list with assay metadata including dimensions and available layers
#' @keywords internal
#' @export
#' @importFrom Seurat VariableFeatures
#' @importFrom SeuratObject Layers
get_assay_info <- function(obj, assay_name) {
  assay <- obj@assays[[assay_name]]

  # Use Layers() to check existence without loading data - much faster for large objects
  available_layers <- SeuratObject::Layers(assay)

  n_variable_features <- length(SeuratObject::VariableFeatures(assay))
  # Safely check for feature metadata - slot may not exist for all assay types (e.g., SCTAssay)
  n_feature_meta_cols <- tryCatch({
    feature_meta <- assay@meta.data
    if (!is.null(feature_meta)) ncol(feature_meta) else 0
  }, error = function(e) 0)

  info <- list(
    name = assay_name,
    features = nrow(assay),
    cells = ncol(assay),
    has_counts = "counts" %in% available_layers,
    has_data = "data" %in% available_layers,
    has_scale_data = "scale.data" %in% available_layers,
    variable_features = n_variable_features,
    has_variable_features = n_variable_features > 0,
    has_feature_metadata = n_feature_meta_cols > 0
  )

  return(info)
}

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
  if (nrow(matrix) > max_rows) {
    row_idx <- 1:max_rows
  } else {
    row_idx <- 1:nrow(matrix)
  }

  if (ncol(matrix) > max_cols) {
    col_idx <- 1:max_cols
  } else {
    col_idx <- 1:ncol(matrix)
  }

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
  if (inherits(matrix, "dgCMatrix") || inherits(matrix, "sparseMatrix")) {
    total_elements <- as.numeric(nrow(matrix)) * as.numeric(ncol(matrix))
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
  } else {
    total_elements <- length(matrix)
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

#' Extract Embeddings
#'
#' Extract dimensional reduction embeddings from a Seurat object.
#'
#' @param obj A Seurat object
#' @param reduction_name Name of the reduction
#' @param dims Vector of two dimensions to extract
#' @return A data frame with cell names and embedding coordinates
#' @keywords internal
#' @export
#' @importFrom SeuratObject Embeddings
extract_embeddings <- function(obj, reduction_name, dims = c(1, 2)) {
  embeddings <- Embeddings(obj, reduction = reduction_name)

  if (max(dims) > ncol(embeddings)) {
    stop("Requested dimensions exceed available dimensions in reduction")
  }

  data.frame(
    Cell = rownames(embeddings),
    Dim1 = embeddings[, dims[1]],
    Dim2 = embeddings[, dims[2]]
  )
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
  errors <- c()

  if (!inherits(obj, "Seurat")) {
    errors <- c(errors, "Object is not a Seurat object")
  }

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

  # get the dimension string for a slot in an assay object
  get_layer_dim <- function(assay_obj, slot) {
    slotnames <- slotNames(assay_obj)
    layer <- NULL
    if ("layers" %in% slotnames) {   # v5 layers slot
      layer <- assay_obj@layers[[slot]]
    }
    else if (slot %in% slotnames) {  # v4 direct slot
      layer <- slot(assay_obj, slot)
    }

    if (is.null(layer)) {
      dim_string <- "0x0"
    } else {
      dims <- dim(layer)
      dim_string <- paste0(dims[1], "x", dims[2])
    }

    return(dim_string)
  }

  slotinfo <- list()
  slots <- c("counts", "data", "scale.data")

  for (assay in assays) {
    assay_obj <- seurat@assays[[assay]]

    defaultassay <- ifelse(default_assay == assay, "YES", "")

    dims <- sapply(slots, function(slot) get_layer_dim(assay_obj, slot))

    hvgs <- length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] <- c(defaultassay, dims, hvgs)
  }

  output$assays_table <- data.frame(do.call(rbind, slotinfo))
  colnames(output$assays_table) <- c("default", slots, "HVGs")

  return(output)
}
