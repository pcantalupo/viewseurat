assay_panel_ui <- function(assay_name, obj) {
  assay_info <- get_assay_info(obj, assay_name)
  
  fluidRow(
    column(12,
      box(
        title = "Assay Information",
        status = "info",
        solidHeader = TRUE,
        width = 12,
        collapsible = TRUE,
        fluidRow(
          column(4, 
            h4("Features (Genes):"), 
            p(strong(format(assay_info$features, big.mark = ",")))
          ),
          column(4, 
            h4("Cells:"), 
            p(strong(format(assay_info$cells, big.mark = ",")))
          ),
          column(4, 
            h4("Variable Features:"), 
            p(strong(format(assay_info$variable_features, big.mark = ",")))
          )
        ),
        fluidRow(
          column(4,
            h5("Has Counts:", if(assay_info$has_counts) "✓" else "✗")
          ),
          column(4,
            h5("Has Data:", if(assay_info$has_data) "✓" else "✗")
          ),
          column(4,
            h5("Has Scale.Data:", if(assay_info$has_scale_data) "✓" else "✗")
          )
        )
      )
    ),
    column(12,
      tabBox(
        id = paste0("assay_", assay_name, "_tabs"),
        width = 12,
        tabPanel("Counts",
          fluidRow(
            column(12,
              if (assay_info$has_counts) {
                list(
                  h4("Counts Matrix"),
                  verbatimTextOutput(paste0(assay_name, "_counts_info")),
                  DTOutput(paste0(assay_name, "_counts_table"))
                )
              } else {
                p("No counts matrix available for this assay.")
              }
            )
          )
        ),
        tabPanel("Data",
          fluidRow(
            column(12,
              if (assay_info$has_data) {
                list(
                  h4("Data (Normalized) Matrix"),
                  verbatimTextOutput(paste0(assay_name, "_data_info")),
                  DTOutput(paste0(assay_name, "_data_table"))
                )
              } else {
                p("No data matrix available for this assay.")
              }
            )
          )
        ),
        tabPanel("Scale.Data",
          fluidRow(
            column(12,
              if (assay_info$has_scale_data) {
                list(
                  h4("Scale.Data Matrix"),
                  verbatimTextOutput(paste0(assay_name, "_scale_info")),
                  DTOutput(paste0(assay_name, "_scale_table"))
                )
              } else {
                p("No scale.data matrix available for this assay.")
              }
            )
          )
        ),
        tabPanel("Variable Features",
          fluidRow(
            column(12,
              h4("Variable Features"),
              DTOutput(paste0(assay_name, "_variable_features")),
              br(),
              plotOutput(paste0(assay_name, "_heatmap"), height = "600px")
            )
          )
        ),
        tabPanel("Feature Metadata",
          fluidRow(
            column(12,
              h4("Feature-level Metadata"),
              DTOutput(paste0(assay_name, "_feature_meta"))
            )
          )
        )
      )
    )
  )
}

assay_panel_server <- function(assay_name, obj, config, output) {
  assay_info <- get_assay_info(obj, assay_name)
  
  if (assay_info$has_counts) {
    output[[paste0(assay_name, "_counts_info")]] <- renderPrint({
      counts_matrix <- get_assay_data_safe(obj, assay_name, "counts")
      if (!is.null(counts_matrix)) {
        sparsity <- get_sparsity_info(counts_matrix)
        cat("Matrix dimensions:", nrow(counts_matrix), "x", ncol(counts_matrix), "\n")
        cat("Sparsity:", sparsity$sparsity_percent, "%\n")
        cat("Memory:", round(sparsity$memory_mb, 2), "MB\n")
        cat("Showing sample of", config$default_matrix_rows, "features and", 
            config$default_matrix_cols, "cells\n")
      }
    })
    
    output[[paste0(assay_name, "_counts_table")]] <- renderDT({
      counts_matrix <- get_assay_data_safe(obj, assay_name, "counts")
      if (!is.null(counts_matrix)) {
        sample_matrix <- get_matrix_sample(
          counts_matrix, 
          config$default_matrix_rows, 
          config$default_matrix_cols
        )
        datatable(
          as.data.frame(sample_matrix),
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            scrollY = "400px"
          )
        )
      }
    })
  }
  
  if (assay_info$has_data) {
    output[[paste0(assay_name, "_data_info")]] <- renderPrint({
      data_matrix <- get_assay_data_safe(obj, assay_name, "data")
      if (!is.null(data_matrix)) {
        sparsity <- get_sparsity_info(data_matrix)
        cat("Matrix dimensions:", nrow(data_matrix), "x", ncol(data_matrix), "\n")
        cat("Sparsity:", sparsity$sparsity_percent, "%\n")
        cat("Memory:", round(sparsity$memory_mb, 2), "MB\n")
        cat("Showing sample of", config$default_matrix_rows, "features and", 
            config$default_matrix_cols, "cells\n")
      }
    })
    
    output[[paste0(assay_name, "_data_table")]] <- renderDT({
      data_matrix <- get_assay_data_safe(obj, assay_name, "data")
      if (!is.null(data_matrix)) {
        sample_matrix <- get_matrix_sample(
          data_matrix, 
          config$default_matrix_rows, 
          config$default_matrix_cols
        )
        datatable(
          as.data.frame(sample_matrix),
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            scrollY = "400px"
          )
        )
      }
    })
  }
  
  if (assay_info$has_scale_data) {
    output[[paste0(assay_name, "_scale_info")]] <- renderPrint({
      scale_matrix <- get_assay_data_safe(obj, assay_name, "scale.data")
      if (!is.null(scale_matrix)) {
        cat("Matrix dimensions:", nrow(scale_matrix), "x", ncol(scale_matrix), "\n")
        cat("Memory:", round(as.numeric(object.size(scale_matrix)) / 1024^2, 2), "MB\n")
        cat("Showing sample of", config$default_matrix_rows, "features and", 
            config$default_matrix_cols, "cells\n")
      }
    })
    
    output[[paste0(assay_name, "_scale_table")]] <- renderDT({
      scale_matrix <- get_assay_data_safe(obj, assay_name, "scale.data")
      if (!is.null(scale_matrix)) {
        sample_matrix <- get_matrix_sample(
          scale_matrix, 
          config$default_matrix_rows, 
          config$default_matrix_cols
        )
        datatable(
          as.data.frame(sample_matrix),
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            scrollY = "400px"
          )
        )
      }
    })
  }
  
  output[[paste0(assay_name, "_variable_features")]] <- renderDT({
    var_features <- get_top_variable_features(obj, assay_name, n = config$max_features_display)
    if (!is.null(var_features) && length(var_features) > 0) {
      datatable(
        data.frame(Feature = var_features, Rank = 1:length(var_features)),
        options = list(pageLength = 25)
      )
    }
  })
  
  output[[paste0(assay_name, "_heatmap")]] <- renderPlot({
    tryCatch({
      plot_feature_heatmap(obj, assay_name, n_cells = 100, config = config)
    }, error = function(e) {
      NULL
    })
  })
  
  output[[paste0(assay_name, "_feature_meta")]] <- renderDT({
    tryCatch({
      assay <- obj@assays[[assay_name]]
      if (nrow(assay@meta.data) > 0) {
        datatable(
          assay@meta.data,
          options = list(
            pageLength = config$rows_per_page,
            scrollX = TRUE
          )
        )
      } else {
        datatable(
          data.frame(Message = "No feature metadata available"),
          options = list(dom = 't')
        )
      }
    }, error = function(e) {
      datatable(
        data.frame(Message = "Unable to load feature metadata"),
        options = list(dom = 't')
      )
    })
  })
}
