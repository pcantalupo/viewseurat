#' Plot Reduction Interactive
#'
#' Create an interactive plotly scatter plot for dimensional reduction.
#'
#' @param obj A Seurat object
#' @param reduction_name Name of the reduction to plot
#' @param color_by Metadata column to color by, or "None"
#' @param dims Vector of two dimensions to plot
#' @return A plotly object
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
plot_reduction_interactive <- function(obj, reduction_name, color_by = "None",
                                       dims = c(1, 2)) {
  embeddings <- extract_embeddings(obj, reduction_name, dims)

  # Sample before joining to avoid copying all metadata for all cells
  if (nrow(embeddings) > 10000) {
    set.seed(42)
    sample_idx <- sample(nrow(embeddings), 10000)
    embeddings <- embeddings[sample_idx, , drop = FALSE]
  }

  # Only join the single column needed for coloring
  if (color_by != "None" && color_by %in% colnames(obj@meta.data)) {
    embeddings[[color_by]] <- obj@meta.data[embeddings$Cell, color_by]
  }
  plot_data <- embeddings

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dim1, y = Dim2, text = Cell)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0(reduction_name, "_", dims[1]),
      y = paste0(reduction_name, "_", dims[2]),
      title = paste(toupper(reduction_name), "Plot")
    )

  if (color_by != "None" && color_by %in% colnames(obj@meta.data)) {
    if (is.numeric(obj@meta.data[[color_by]])) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes_string(color = color_by),
        size = 1.0,
        alpha = 0.8
      ) +
        ggplot2::scale_color_viridis_c(option = "viridis")
    } else {
      p <- p + ggplot2::geom_point(
        ggplot2::aes_string(color = color_by),
        size = 1.0,
        alpha = 0.8
      )
    }
  } else {
    p <- p + ggplot2::geom_point(
      size = 1.0,
      alpha = 0.8,
      color = "steelblue"
    )
  }

  plotly::ggplotly(p, tooltip = c("text", "x", "y", if(color_by != "None") color_by)) %>%
    plotly::layout(hovermode = 'closest')
}

#' Plot Reduction Static
#'
#' Create a static ggplot scatter plot for dimensional reduction.
#' Faster than interactive version for large datasets.
#'
#' @param obj A Seurat object
#' @param reduction_name Name of the reduction to plot
#' @param color_by Metadata column to color by, or "None"
#' @param dims Vector of two dimensions to plot
#' @return A ggplot object
#' @keywords internal
#' @export
plot_reduction_static <- function(obj, reduction_name, color_by = "None",
                                  dims = c(1, 2)) {
  embeddings <- extract_embeddings(obj, reduction_name, dims)

  # Sample before joining to avoid copying all metadata for all cells
  if (nrow(embeddings) > 10000) {
    set.seed(42)
    sample_idx <- sample(nrow(embeddings), 10000)
    embeddings <- embeddings[sample_idx, , drop = FALSE]
  }

  # Only join the single column needed for coloring
  if (color_by != "None" && color_by %in% colnames(obj@meta.data)) {
    embeddings[[color_by]] <- obj@meta.data[embeddings$Cell, color_by]
  }
  plot_data <- embeddings

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dim1, y = Dim2)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0(reduction_name, "_", dims[1]),
      y = paste0(reduction_name, "_", dims[2]),
      title = paste(toupper(reduction_name), "Plot")
    )

  if (color_by != "None" && color_by %in% colnames(obj@meta.data)) {
    if (is.numeric(obj@meta.data[[color_by]])) {
      p <- p + ggplot2::geom_point(
        ggplot2::aes_string(color = color_by),
        size = 1.0,
        alpha = 0.8
      ) +
        ggplot2::scale_color_viridis_c(option = "viridis")
    } else {
      p <- p + ggplot2::geom_point(
        ggplot2::aes_string(color = color_by),
        size = 1.0,
        alpha = 0.8
      )
    }
  } else {
    p <- p + ggplot2::geom_point(
      size = 1.0,
      alpha = 0.8,
      color = "steelblue"
    )
  }

  return(p)
}

#' Plot Metadata Distribution
#'
#' Create a histogram or bar plot showing distribution of a metadata column.
#'
#' @param obj A Seurat object
#' @param column Metadata column name
#' @return A ggplot object
#' @keywords internal
#' @export
plot_metadata_distribution <- function(obj, column) {
  data <- obj@meta.data[[column]]

  if (is.numeric(data)) {
    p <- ggplot2::ggplot(obj@meta.data, ggplot2::aes_string(x = column)) +
      ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Distribution of", column),
        x = column,
        y = "Count"
      ) +
      ggplot2::theme(
        axis.text = ggplot2::element_text(size = 14, color = "black"),
        axis.title = ggplot2::element_text(size = 16, color = "black", face = "bold"),
        plot.title = ggplot2::element_text(size = 18, color = "black", face = "bold")
      )
  } else {
    freq_table <- as.data.frame(table(data))
    colnames(freq_table) <- c("Category", "Count")

    p <- ggplot2::ggplot(freq_table, ggplot2::aes(x = Category, y = Count)) +
      ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("Distribution of", column),
        x = column,
        y = "Count"
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 14, color = "black"),
        axis.text.y = ggplot2::element_text(size = 14, color = "black"),
        axis.title = ggplot2::element_text(size = 16, color = "black", face = "bold"),
        plot.title = ggplot2::element_text(size = 18, color = "black", face = "bold")
      )
  }

  return(p)
}


#' Plot QC Metrics
#'
#' Create histograms of common QC metrics.
#'
#' @param obj A Seurat object
#' @return A list of ggplot objects or NULL
#' @keywords internal
#' @export
plot_qc_metrics <- function(obj) {
  qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  available_cols <- intersect(qc_cols, colnames(obj@meta.data))

  if (length(available_cols) == 0) {
    return(NULL)
  }

  plots <- list()

  for (col in available_cols) {
    p <- ggplot2::ggplot(obj@meta.data, ggplot2::aes_string(x = col)) +
      ggplot2::geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = median(obj@meta.data[[col]])),
                 color = "red", linetype = "dashed") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = col, x = col, y = "Count")

    plots[[col]] <- p
  }

  return(plots)
}

#' Create Matrix Heatmap
#'
#' Create an interactive heatmap of matrix data.
#'
#' @param matrix_data A matrix to plot
#' @param title Plot title
#' @return A plotly object
#' @keywords internal
#' @export
create_matrix_heatmap <- function(matrix_data, title = "Matrix Heatmap") {
  if (inherits(matrix_data, "dgCMatrix") || inherits(matrix_data, "sparseMatrix")) {
    matrix_data <- as.matrix(matrix_data)
  }

  plotly::plot_ly(
    z = matrix_data,
    x = colnames(matrix_data),
    y = rownames(matrix_data),
    type = "heatmap",
    colorscale = "Viridis"
  ) %>%
    plotly::layout(
      title = title,
      xaxis = list(title = "Cells", tickangle = -45),
      yaxis = list(title = "Features")
    )
}
