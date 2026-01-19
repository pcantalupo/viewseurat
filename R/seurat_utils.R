#' Load Configuration
#'
#' Load configuration from config.yaml, merging with default values.
#' Searches for config.yaml in the current directory first, then in the
#' package's inst/app directory.
#'
#' @return A list of configuration values
#' @export
load_config <- function() {
  default_config <- list(
    max_cells_display = 10000,
    max_features_display = 100,
    matrix_sample_size = 1000,
    default_reduction = "umap",
    color_palette = "viridis",
    plot_theme = "minimal",
    point_size = 1.0,
    point_alpha = 0.8,
    show_legend = TRUE,
    show_sparse_info = TRUE,
    default_matrix_rows = 10,
    default_matrix_cols = 20,
    enable_caching = TRUE,
    parallel_processing = FALSE,
    max_upload_size_mb = 10240,
    rows_per_page = 25
  )

  # Try to find config.yaml in various locations
  config_paths <- c(
    "config.yaml",  # Current directory (for standalone use)
    system.file("app", "config.yaml", package = "viewseurat")  # Package installation
  )

  config_file <- NULL
  for (path in config_paths) {
    if (path != "" && file.exists(path)) {
      config_file <- path
      break
    }
  }

  if (!is.null(config_file)) {
    tryCatch({
      user_config <- yaml::read_yaml(config_file)
      for (name in names(user_config)) {
        default_config[[name]] <- user_config[[name]]
      }
      message("Loaded custom configuration from ", config_file)
    }, error = function(e) {
      warning("Could not load config.yaml, using defaults: ", e$message)
    })
  }

  return(default_config)
}

#' Get Assay Information
#'
#' Extract information about a specific assay from a Seurat object.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay to get information for
#' @return A list with assay metadata including dimensions and available layers
#' @keywords internal
#' @export
#' @importFrom Seurat GetAssayData VariableFeatures
get_assay_info <- function(obj, assay_name) {
  assay <- obj@assays[[assay_name]]

  info <- list(
    name = assay_name,
    features = nrow(assay),
    cells = ncol(assay),
    has_counts = length(GetAssayData(assay, layer = "counts")) > 0,
    has_data = length(GetAssayData(assay, layer = "data")) > 0,
    has_scale_data = length(GetAssayData(assay, layer = "scale.data")) > 0,
    variable_features = length(VariableFeatures(assay))
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

#' Get Assay Data Safely
#'
#' Wrapper around GetAssayData with error handling.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay
#' @param layer Name of the layer (counts, data, scale.data)
#' @return Matrix data or NULL if not available
#' @keywords internal
#' @export
#' @importFrom Seurat GetAssayData DefaultAssay
get_assay_data_safe <- function(obj, assay_name, layer) {
  tryCatch({
    DefaultAssay(obj) <- assay_name
    GetAssayData(obj, layer = layer)
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
