plot_reduction_interactive <- function(obj, reduction_name, color_by = "None", 
                                       dims = c(1, 2), config) {
  embeddings <- extract_embeddings(obj, reduction_name, dims)
  
  plot_data <- cbind(embeddings, obj@meta.data)
  
  p <- ggplot(plot_data, aes(x = Dim1, y = Dim2, text = Cell)) +
    theme_minimal() +
    labs(
      x = paste0(reduction_name, "_", dims[1]),
      y = paste0(reduction_name, "_", dims[2]),
      title = paste(toupper(reduction_name), "Plot")
    )
  
  if (color_by != "None" && color_by %in% colnames(obj@meta.data)) {
    if (is.numeric(obj@meta.data[[color_by]])) {
      p <- p + geom_point(
        aes_string(color = color_by),
        size = config$point_size,
        alpha = config$point_alpha
      ) +
        scale_color_viridis_c(option = config$color_palette)
    } else {
      p <- p + geom_point(
        aes_string(color = color_by),
        size = config$point_size,
        alpha = config$point_alpha
      )
    }
  } else {
    p <- p + geom_point(
      size = config$point_size,
      alpha = config$point_alpha,
      color = "steelblue"
    )
  }
  
  ggplotly(p, tooltip = c("text", "x", "y", if(color_by != "None") color_by)) %>%
    layout(hovermode = 'closest')
}

plot_feature_expression <- function(obj, reduction_name, feature, dims = c(1, 2), config) {
  embeddings <- extract_embeddings(obj, reduction_name, dims)
  
  if (!feature %in% rownames(obj)) {
    stop("Feature not found in object")
  }
  
  expression <- GetAssayData(obj, layer = "data")[feature, ]
  plot_data <- cbind(embeddings, Expression = expression)
  
  p <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = Expression)) +
    geom_point(size = config$point_size, alpha = config$point_alpha) +
    scale_color_viridis_c(option = config$color_palette) +
    theme_minimal() +
    labs(
      x = paste0(reduction_name, "_", dims[1]),
      y = paste0(reduction_name, "_", dims[2]),
      title = paste("Expression of", feature)
    )
  
  return(p)
}

plot_metadata_distribution <- function(obj, column, config) {
  data <- obj@meta.data[[column]]
  
  if (is.numeric(data)) {
    p <- ggplot(obj@meta.data, aes_string(x = column)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      labs(
        title = paste("Distribution of", column),
        x = column,
        y = "Count"
      )
  } else {
    freq_table <- as.data.frame(table(data))
    colnames(freq_table) <- c("Category", "Count")
    
    p <- ggplot(freq_table, aes(x = Category, y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      theme_minimal() +
      labs(
        title = paste("Distribution of", column),
        x = column,
        y = "Count"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}

plot_feature_heatmap <- function(obj, assay_name, features = NULL, n_cells = 100, config) {
  DefaultAssay(obj) <- assay_name
  
  if (is.null(features)) {
    features <- get_top_variable_features(obj, assay_name, n = 20)
  }
  
  if (is.null(features) || length(features) == 0) {
    return(NULL)
  }
  
  if (ncol(obj) > n_cells) {
    cells_to_plot <- sample(colnames(obj), n_cells)
  } else {
    cells_to_plot <- colnames(obj)
  }
  
  data_matrix <- GetAssayData(obj, layer = "scale.data")
  
  if (nrow(data_matrix) == 0) {
    data_matrix <- GetAssayData(obj, layer = "data")
  }
  
  features <- intersect(features, rownames(data_matrix))
  if (length(features) == 0) {
    return(NULL)
  }
  
  plot_data <- as.matrix(data_matrix[features, cells_to_plot, drop = FALSE])
  
  p <- DoHeatmap(
    obj,
    features = features,
    cells = cells_to_plot,
    size = 3
  ) +
    theme(axis.text.y = element_text(size = 8))
  
  return(p)
}

plot_qc_metrics <- function(obj, config) {
  qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  available_cols <- intersect(qc_cols, colnames(obj@meta.data))
  
  if (length(available_cols) == 0) {
    return(NULL)
  }
  
  plots <- list()
  
  for (col in available_cols) {
    p <- ggplot(obj@meta.data, aes_string(x = col)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      geom_vline(aes(xintercept = median(obj@meta.data[[col]])), 
                 color = "red", linetype = "dashed") +
      theme_minimal() +
      labs(title = col, x = col, y = "Count")
    
    plots[[col]] <- p
  }
  
  return(plots)
}

create_matrix_heatmap <- function(matrix_data, title = "Matrix Heatmap") {
  if (inherits(matrix_data, "dgCMatrix") || inherits(matrix_data, "sparseMatrix")) {
    matrix_data <- as.matrix(matrix_data)
  }
  
  plot_ly(
    z = matrix_data,
    x = colnames(matrix_data),
    y = rownames(matrix_data),
    type = "heatmap",
    colorscale = "Viridis"
  ) %>%
    layout(
      title = title,
      xaxis = list(title = "Cells", tickangle = -45),
      yaxis = list(title = "Features")
    )
}
