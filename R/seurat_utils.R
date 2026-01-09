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
    default_matrix_rows = 50,
    default_matrix_cols = 20,
    enable_caching = TRUE,
    parallel_processing = FALSE,
    max_upload_size_mb = 10240,
    rows_per_page = 25
  )
  
  config_file <- "config.yaml"
  if (file.exists(config_file)) {
    tryCatch({
      library(yaml)
      user_config <- yaml::read_yaml(config_file)
      for (name in names(user_config)) {
        default_config[[name]] <- user_config[[name]]
      }
      message("Loaded custom configuration from config.yaml")
    }, error = function(e) {
      warning("Could not load config.yaml, using defaults: ", e$message)
    })
  }
  
  return(default_config)
}

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

get_matrix_sample <- function(matrix, max_rows = 50, max_cols = 20) {
  if (nrow(matrix) > max_rows) {
    row_idx <- sample(nrow(matrix), max_rows)
  } else {
    row_idx <- 1:nrow(matrix)
  }
  
  if (ncol(matrix) > max_cols) {
    col_idx <- sample(ncol(matrix), max_cols)
  } else {
    col_idx <- 1:ncol(matrix)
  }
  
  if (inherits(matrix, "dgCMatrix") || inherits(matrix, "sparseMatrix")) {
    return(as.matrix(matrix[row_idx, col_idx, drop = FALSE]))
  } else {
    return(matrix[row_idx, col_idx, drop = FALSE])
  }
}

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

get_top_variable_features <- function(obj, assay_name, n = 20) {
  DefaultAssay(obj) <- assay_name
  var_features <- VariableFeatures(obj)
  
  if (length(var_features) == 0) {
    return(NULL)
  }
  
  head(var_features, n)
}

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

get_assay_data_safe <- function(obj, assay_name, layer) {
  tryCatch({
    DefaultAssay(obj) <- assay_name
    GetAssayData(obj, layer = layer)
  }, error = function(e) {
    NULL
  })
}
