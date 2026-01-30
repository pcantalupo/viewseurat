#' Assay Panel UI
#'
#' Generate the UI for a single assay panel.
#'
#' @param assay_name Name of the assay
#' @param obj A Seurat object
#' @return A shiny fluidRow UI element
#' @keywords internal
#' @export
assay_panel_ui <- function(assay_name, obj) {
  assay_info <- get_assay_info(obj, assay_name)

  shiny::fluidRow(
    shiny::column(12,
      shinydashboard::box(
        title = "Assay Information",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        shiny::fluidRow(
          shiny::column(4,
            shiny::h4("Features (Genes):"),
            shiny::p(shiny::strong(format(assay_info$features, big.mark = ",")))
          ),
          shiny::column(4,
            shiny::h4("Cells:"),
            shiny::p(shiny::strong(format(assay_info$cells, big.mark = ",")))
          ),
          shiny::column(4,
            shiny::h4("Variable Features:"),
            shiny::p(shiny::strong(format(assay_info$variable_features, big.mark = ",")))
          )
        ),
      )
    ),
    shiny::column(12,
      shinydashboard::tabBox(
        id = paste0("assay_", assay_name, "_tabs"),
        width = 12,
        shiny::tabPanel(
          if (assay_info$has_counts) "Counts \u2713" else "Counts",
          shiny::fluidRow(
            shiny::column(12,
              if (assay_info$has_counts) {
                list(
                  shiny::h4("Counts Matrix"),
                  shiny::verbatimTextOutput(paste0(assay_name, "_counts_info")),
                  shinycssloaders::withSpinner(DT::DTOutput(paste0(assay_name, "_counts_table")))
                )
              } else {
                shiny::p("No counts matrix available for this assay.")
              }
            )
          )
        ),
        shiny::tabPanel(
          if (assay_info$has_data) "Data \u2713" else "Data",
          shiny::fluidRow(
            shiny::column(12,
              if (assay_info$has_data) {
                list(
                  shiny::h4("Data (Normalized) Matrix"),
                  shiny::verbatimTextOutput(paste0(assay_name, "_data_info")),
                  shinycssloaders::withSpinner(DT::DTOutput(paste0(assay_name, "_data_table")))
                )
              } else {
                shiny::p("No data matrix available for this assay.")
              }
            )
          )
        ),
        shiny::tabPanel(
          if (assay_info$has_scale_data) "Scale.Data \u2713" else "Scale.Data",
          shiny::fluidRow(
            shiny::column(12,
              if (assay_info$has_scale_data) {
                list(
                  shiny::h4("Scale.Data Matrix"),
                  shiny::verbatimTextOutput(paste0(assay_name, "_scale_info")),
                  shinycssloaders::withSpinner(DT::DTOutput(paste0(assay_name, "_scale_table")))
                )
              } else {
                shiny::p("No scale.data matrix available for this assay.")
              }
            )
          )
        ),
        shiny::tabPanel(
          if (assay_info$has_variable_features) "Variable Features \u2713" else "Variable Features",
          shiny::fluidRow(
            shiny::column(12,
              if (assay_info$has_variable_features) {
                list(
                  shiny::h4("Variable Features"),
                  DT::DTOutput(paste0(assay_name, "_variable_features"))
                )
              } else {
                shiny::p("No variable features for this assay.")
              }
            )
          )
        ),
        shiny::tabPanel(
          if (assay_info$has_feature_metadata) "Feature Metadata \u2713" else "Feature Metadata",
          shiny::fluidRow(
            shiny::column(12,
              if (assay_info$has_feature_metadata) {
                list(
                  shiny::h4("Feature-level Metadata"),
                  DT::DTOutput(paste0(assay_name, "_feature_meta"))
                )
              } else {
                shiny::p("No feature metadata for this assay.")
              }
            )
          )
        )
      )
    )
  )
}

#' Assay Panel Server
#'
#' Server logic for rendering assay data.
#'
#' @param assay_name Name of the assay
#' @param obj A Seurat object
#' @param config Configuration list
#' @param output Shiny output object
#' @return NULL (called for side effects)
#' @keywords internal
#' @export
assay_panel_server <- function(assay_name, obj, config, output) {
  assay_info <- get_assay_info(obj, assay_name)

  if (assay_info$has_counts) {
    output[[paste0(assay_name, "_counts_info")]] <- shiny::renderPrint({
      counts_matrix <- get_assay_data_safe(obj, assay_name, "counts")
      if (!is.null(counts_matrix)) {
        sparsity <- get_sparsity_info(counts_matrix)
        cat("Matrix dimensions:", nrow(counts_matrix), "x", ncol(counts_matrix), "\n")
        cat("Sparsity:", sparsity$sparsity_percent, "%\n")
        cat("Memory:", round(sparsity$memory_mb, 2), "MB\n")
        cat("Showing first", min(nrow(counts_matrix), config$default_matrix_rows),
            "of", nrow(counts_matrix), "features and first",
            min(ncol(counts_matrix), config$default_matrix_cols), "cells\n")
      }
    })

    output[[paste0(assay_name, "_counts_table")]] <- DT::renderDT({
      counts_matrix <- get_assay_data_safe(obj, assay_name, "counts")
      if (!is.null(counts_matrix)) {
        # Limit rows to avoid expensive sparse-to-dense conversion of all genes
        max_rows <- min(nrow(counts_matrix), config$default_matrix_rows)
        num_cols <- min(ncol(counts_matrix), config$default_matrix_cols)
        row_idx <- 1:max_rows
        col_idx <- 1:num_cols

        if (inherits(counts_matrix, "dgCMatrix") || inherits(counts_matrix, "sparseMatrix")) {
          sample_matrix <- as.matrix(counts_matrix[row_idx, col_idx, drop = FALSE])
        } else {
          sample_matrix <- counts_matrix[row_idx, col_idx, drop = FALSE]
        }

        DT::datatable(
          as.data.frame(sample_matrix),
          options = list(
            pageLength = 10,
            lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
            scrollX = TRUE,
            scrollY = "400px"
          )
        )
      }
    }, server = TRUE)
  }

  if (assay_info$has_data) {
    output[[paste0(assay_name, "_data_info")]] <- shiny::renderPrint({
      data_matrix <- get_assay_data_safe(obj, assay_name, "data")
      if (!is.null(data_matrix)) {
        sparsity <- get_sparsity_info(data_matrix)
        cat("Matrix dimensions:", nrow(data_matrix), "x", ncol(data_matrix), "\n")
        cat("Sparsity:", sparsity$sparsity_percent, "%\n")
        cat("Memory:", round(sparsity$memory_mb, 2), "MB\n")
        cat("Showing first", min(nrow(data_matrix), config$default_matrix_rows),
            "of", nrow(data_matrix), "features and first",
            min(ncol(data_matrix), config$default_matrix_cols), "cells\n")
      }
    })

    output[[paste0(assay_name, "_data_table")]] <- DT::renderDT({
      data_matrix <- get_assay_data_safe(obj, assay_name, "data")
      if (!is.null(data_matrix)) {
        # Limit rows to avoid expensive sparse-to-dense conversion of all genes
        max_rows <- min(nrow(data_matrix), config$default_matrix_rows)
        num_cols <- min(ncol(data_matrix), config$default_matrix_cols)
        row_idx <- 1:max_rows
        col_idx <- 1:num_cols

        if (inherits(data_matrix, "dgCMatrix") || inherits(data_matrix, "sparseMatrix")) {
          sample_matrix <- as.matrix(data_matrix[row_idx, col_idx, drop = FALSE])
        } else {
          sample_matrix <- data_matrix[row_idx, col_idx, drop = FALSE]
        }

        DT::datatable(
          as.data.frame(sample_matrix),
          options = list(
            pageLength = 10,
            lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
            scrollX = TRUE,
            scrollY = "400px"
          )
        )
      }
    }, server = TRUE)
  }

  if (assay_info$has_scale_data) {
    output[[paste0(assay_name, "_scale_info")]] <- shiny::renderPrint({
      scale_matrix <- get_assay_data_safe(obj, assay_name, "scale.data")
      if (!is.null(scale_matrix)) {
        cat("Matrix dimensions:", nrow(scale_matrix), "x", ncol(scale_matrix), "\n")
        cat("Memory:", round(as.numeric(object.size(scale_matrix)) / 1024^2, 2), "MB\n")
        cat("Showing first", min(nrow(scale_matrix), config$default_matrix_rows),
            "of", nrow(scale_matrix), "features and first",
            min(ncol(scale_matrix), config$default_matrix_cols), "cells\n")
      }
    })

    output[[paste0(assay_name, "_scale_table")]] <- DT::renderDT({
      scale_matrix <- get_assay_data_safe(obj, assay_name, "scale.data")
      if (!is.null(scale_matrix)) {
        # Limit rows to avoid expensive sparse-to-dense conversion of all genes
        max_rows <- min(nrow(scale_matrix), config$default_matrix_rows)
        num_cols <- min(ncol(scale_matrix), config$default_matrix_cols)
        row_idx <- 1:max_rows
        col_idx <- 1:num_cols

        if (inherits(scale_matrix, "dgCMatrix") || inherits(scale_matrix, "sparseMatrix")) {
          sample_matrix <- as.matrix(scale_matrix[row_idx, col_idx, drop = FALSE])
        } else {
          sample_matrix <- scale_matrix[row_idx, col_idx, drop = FALSE]
        }

        DT::datatable(
          as.data.frame(sample_matrix),
          options = list(
            pageLength = 10,
            lengthMenu = list(c(10, 25, 50, 100, -1), c('10', '25', '50', '100', 'All')),
            scrollX = TRUE,
            scrollY = "400px"
          )
        )
      }
    }, server = TRUE)
  }

  output[[paste0(assay_name, "_variable_features")]] <- DT::renderDT({
    var_features <- get_top_variable_features(obj, assay_name, n = config$max_features_display)
    if (!is.null(var_features) && length(var_features) > 0) {
      var_features_sorted <- sort(var_features)
      DT::datatable(
        data.frame(Feature = var_features_sorted),
        options = list(pageLength = 25)
      )
    }
  })

  output[[paste0(assay_name, "_feature_meta")]] <- DT::renderDT({
    tryCatch({
      assay <- obj@assays[[assay_name]]
      if (nrow(assay@meta.data) > 0) {
        DT::datatable(
          assay@meta.data,
          options = list(
            pageLength = config$rows_per_page,
            scrollX = TRUE
          )
        )
      } else {
        DT::datatable(
          data.frame(Message = "No feature metadata available"),
          options = list(dom = 't')
        )
      }
    }, error = function(e) {
      DT::datatable(
        data.frame(Message = "No feature metadata for this assay"),
        options = list(dom = 't')
      )
    })
  })
}
