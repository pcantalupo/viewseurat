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
#' @param output Shiny output object
#' @return NULL (called for side effects)
#' @keywords internal
#' @export
assay_panel_server <- function(assay_name, obj, output) {
  assay_info <- get_assay_info(obj, assay_name)


  # Constants for preview limits

max_preview_features <- 50
  max_preview_cells <- 20

  # Cache layer data with reactive() - fetch once with subsetting, reuse for info + table
  # LayerData() with features/cells params subsets at retrieval time,
 # which is much faster than loading the full matrix then subsetting
  counts_data <- shiny::reactive({
    if (assay_info$has_counts) {
      get_assay_data_safe(obj, assay_name, "counts",
                          max_features = max_preview_features,
                          max_cells = max_preview_cells)
    }
  })

  data_data <- shiny::reactive({
    if (assay_info$has_data) {
      get_assay_data_safe(obj, assay_name, "data",
                          max_features = max_preview_features,
                          max_cells = max_preview_cells)
    }
  })

  scale_data <- shiny::reactive({
    if (assay_info$has_scale_data) {
      get_assay_data_safe(obj, assay_name, "scale.data",
                          max_features = max_preview_features,
                          max_cells = max_preview_cells)
    }
  })

  # Helper to convert sparse matrix to dense for display
  to_dense <- function(mat) {
    if (inherits(mat, "dgCMatrix") || inherits(mat, "sparseMatrix")) {
      as.matrix(mat)
    } else {
      mat
    }
  }

  # Helper to create standard DataTable options
  dt_options <- function() {
    list(
      pageLength = 10,
      lengthMenu = list(c(10, 25, 50, -1), c('10', '25', '50', 'All')),
      scrollX = TRUE,
      scrollY = "400px"
    )
  }

  if (assay_info$has_counts) {
    output[[paste0(assay_name, "_counts_info")]] <- shiny::renderPrint({
      d <- get_layer_dim(obj@assays[[assay_name]], "counts")
      full_rows <- d[1]
      full_cols <- d[2]

      cat("Matrix dimensions:", full_rows, "x", full_cols, "\n")
    })

    output[[paste0(assay_name, "_counts_table")]] <- DT::renderDT({
      mat <- counts_data()
      if (!is.null(mat)) {
        DT::datatable(
          as.data.frame(to_dense(mat)),
          options = dt_options()
        )
      }
    }, server = TRUE)
  }

  if (assay_info$has_data) {
    output[[paste0(assay_name, "_data_info")]] <- shiny::renderPrint({
      d <- get_layer_dim(obj@assays[[assay_name]], "data")
      full_rows <- d[1]
      full_cols <- d[2]

      cat("Matrix dimensions:", full_rows, "x", full_cols, "\n")
    })

    output[[paste0(assay_name, "_data_table")]] <- DT::renderDT({
      mat <- data_data()
      if (!is.null(mat)) {
        DT::datatable(
          as.data.frame(to_dense(mat)),
          options = dt_options()
        )
      }
    }, server = TRUE)
  }

  if (assay_info$has_scale_data) {
    output[[paste0(assay_name, "_scale_info")]] <- shiny::renderPrint({
      d <- get_layer_dim(obj@assays[[assay_name]], "scale.data")
      full_rows <- d[1]
      full_cols <- d[2]

      cat("Matrix dimensions:", full_rows, "x", full_cols, "\n")
    })

    output[[paste0(assay_name, "_scale_table")]] <- DT::renderDT({
      mat <- scale_data()
      if (!is.null(mat)) {
        DT::datatable(
          as.data.frame(to_dense(mat)),
          options = dt_options()
        )
      }
    }, server = TRUE)
  }

  output[[paste0(assay_name, "_variable_features")]] <- DT::renderDT({
    var_features <- get_top_variable_features(obj, assay_name, n = 100)
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
            pageLength = 10,
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
