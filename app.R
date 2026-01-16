library(shiny)
library(shinydashboard)
library(Seurat)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(Matrix)
library(qs2)
library(shinyjs)

source("R/seurat_utils.R")
source("R/plot_functions.R")
source("R/ui_modules.R")

config <- load_config()
options(shiny.maxRequestSize = config$max_upload_size_mb * 1024^2)

ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(title = "Seurat Object Viewer"),
  
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
      menuItem("Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Overview", tabName = "overview", icon = icon("info-circle")),
      menuItem("Assays", tabName = "assays", icon = icon("table")),
      menuItem("Reductions", tabName = "reductions", icon = icon("project-diagram")),
      menuItem("Metadata", tabName = "metadata", icon = icon("list")),
      menuItem("Graphs", tabName = "graphs", icon = icon("share-alt")),
      menuItem("Images", tabName = "images", icon = icon("image"))
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
        .box-body { overflow-x: auto; }
        .dataTables_wrapper { overflow-x: auto; }
        .small-box { cursor: pointer; }
        #drop_zone {
          border: 3px dashed #3c8dbc;
          border-radius: 10px;
          padding: 50px;
          text-align: center;
          background-color: #f9f9f9;
          transition: all 0.3s;
          margin-bottom: 20px;
        }
        #drop_zone.dragover {
          border-color: #00a65a;
          background-color: #e8f5e9;
        }
        #drop_zone h3 {
          color: #3c8dbc;
          margin-bottom: 10px;
        }
        #drop_zone p {
          color: #666;
        }
      ")),
      tags$script(HTML("
        $(document).on('shiny:connected', function() {
          var dropZone = document.getElementById('drop_zone');
          
          // Prevent default drag behaviors
          ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(eventName => {
            dropZone.addEventListener(eventName, preventDefaults, false);
            document.body.addEventListener(eventName, preventDefaults, false);
          });
          
          function preventDefaults(e) {
            e.preventDefault();
            e.stopPropagation();
          }
          
          // Highlight drop zone when dragging over it
          ['dragenter', 'dragover'].forEach(eventName => {
            dropZone.addEventListener(eventName, highlight, false);
          });
          
          ['dragleave', 'drop'].forEach(eventName => {
            dropZone.addEventListener(eventName, unhighlight, false);
          });
          
          function highlight(e) {
            dropZone.classList.add('dragover');
          }
          
          function unhighlight(e) {
            dropZone.classList.remove('dragover');
          }
          
          // Handle dropped files
          dropZone.addEventListener('drop', handleDrop, false);
          
          function handleDrop(e) {
            var dt = e.dataTransfer;
            var files = dt.files;
            
            if (files.length > 0) {
              var file = files[0];
              
              // Check if file is .rds or .qs2
              if (file.name.toLowerCase().endsWith('.rds') || file.name.toLowerCase().endsWith('.qs2')) {
                // Update the file input
                var fileInput = document.getElementById('seurat_file');
                var dataTransfer = new DataTransfer();
                dataTransfer.items.add(file);
                fileInput.files = dataTransfer.files;
                
                // Trigger change event
                var event = new Event('change', { bubbles: true });
                fileInput.dispatchEvent(event);
                
                // Update drop zone text
                dropZone.innerHTML = '<h3><i class=\"fa fa-check-circle\"></i> File Ready</h3><p>' + file.name + '</p>';
              } else {
                alert('Please upload an .rds or .qs2 file');
              }
            }
          }
        });
      "))
    ),
    
    tabItems(
      tabItem(tabName = "upload",
        fluidRow(
          box(
            title = "Upload Seurat Object", 
            status = "primary", 
            solidHeader = TRUE,
            width = 12,
            div(
              id = "drop_zone",
              h3(icon("cloud-upload"), "Drag & Drop RDS or QS2 File Here"),
              p("or click below to browse"),
              br()
            ),
            fileInput("seurat_file", 
                      "Choose RDS or QS2 file containing Seurat object",
                      accept = c(".rds", ".RDS", ".qs2", ".QS2")),
            verbatimTextOutput("upload_status")
          )
        )
      ),
      
      tabItem(tabName = "overview",
        uiOutput("overview_ui")
      ),
      
      tabItem(tabName = "assays",
        uiOutput("assays_ui")
      ),
      
      tabItem(tabName = "reductions",
        uiOutput("reductions_ui")
      ),
      
      tabItem(tabName = "metadata",
        uiOutput("metadata_ui")
      ),
      
      tabItem(tabName = "graphs",
        uiOutput("graphs_ui")
      ),
      
      tabItem(tabName = "images",
        uiOutput("images_ui")
      )
    )
  )
)

server <- function(input, output, session) {
  
  seurat_obj <- reactiveVal(NULL)
  
  observeEvent(input$seurat_file, {
    req(input$seurat_file)
    
    tryCatch({
      # Detect file type and load accordingly
      file_ext <- tolower(tools::file_ext(input$seurat_file$name))
      
      if (file_ext == "rds") {
        obj <- readRDS(input$seurat_file$datapath)
      } else if (file_ext == "qs2") {
        obj <- qs2::qs_read(input$seurat_file$datapath)
      } else {
        stop("Unsupported file format. Please use .rds or .qs2 files.")
      }
      
      if (!inherits(obj, "Seurat")) {
        stop("File does not contain a valid Seurat object")
      }
      
      seurat_obj(obj)
      
      output$upload_status <- renderText({
        paste0(
          "✓ Successfully loaded Seurat object\n",
          "File type: ", toupper(file_ext), "\n",
          "Cells: ", ncol(obj), "\n",
          "Assays: ", paste(names(obj@assays), collapse = ", "), "\n",
          "Reductions: ", paste(names(obj@reductions), collapse = ", ")
        )
      })
      
      showNotification("Seurat object loaded successfully!", type = "message")
      updateTabItems(session, "sidebar", "overview")
      
    }, error = function(e) {
      output$upload_status <- renderText({
        paste("✗ Error loading file:", e$message)
      })
      showNotification(paste("Error:", e$message), type = "error", duration = NULL)
    })
  })
  
  output$overview_ui <- renderUI({
    req(seurat_obj())
    obj <- seurat_obj()
    
    fluidRow(
      column(12,
        box(
          title = "Seurat Summary",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          verbatimTextOutput("object_summary")
        )
      ),
      column(12,
        box(
          title = "SeuratInfo",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          verbatimTextOutput("seurat_info_output")
        )
      ),
      column(6,
        valueBox(
          value = ncol(obj),
          subtitle = "Total Cells",
          icon = icon("circle"),
          color = "blue"
        )
      ),
      column(6,
        div(
          onclick = "Shiny.setInputValue('goto_assays', Math.random());",
          style = "cursor: pointer;",
          valueBox(
            value = length(obj@assays),
            subtitle = "Assays",
            icon = icon("table"),
            color = "green"
          )
        )
      ),
      column(6,
        div(
          onclick = "Shiny.setInputValue('goto_reductions', Math.random());",
          style = "cursor: pointer;",
          valueBox(
            value = length(obj@reductions),
            subtitle = "Reductions",
            icon = icon("project-diagram"),
            color = "purple"
          )
        )
      ),
      column(6,
        div(
          onclick = "Shiny.setInputValue('goto_metadata', Math.random());",
          style = "cursor: pointer;",
          valueBox(
            value = ncol(obj@meta.data),
            subtitle = "Metadata Columns",
            icon = icon("list"),
            color = "orange"
          )
        )
      )
    )
  })
  
  output$object_summary <- renderPrint({
    req(seurat_obj())
    print(seurat_obj())
  })
  
  output$seurat_info_output <- renderPrint({
    req(seurat_obj())
    obj <- seurat_obj()
    info <- SeuratInfo(obj)
    
    cat(info$version, "\n")
    cat(info$graphs, "\n")
    cat(info$reductions, "\n")
    cat(info$images, "\n")
    cat("\n", info$ident_label, "\n")
    cat("Idents():\n")
    print(info$idents_table)
    cat("\nAssays:\n")
    print(info$assays_table)
  })
  
  observeEvent(input$goto_assays, {
    updateTabItems(session, "sidebar", "assays")
  })
  
  observeEvent(input$goto_reductions, {
    updateTabItems(session, "sidebar", "reductions")
  })
  
  observeEvent(input$goto_metadata, {
    updateTabItems(session, "sidebar", "metadata")
  })
  
  output$assays_ui <- renderUI({
    req(seurat_obj())
    obj <- seurat_obj()
    assay_names <- names(obj@assays)
    
    do.call(tabBox, c(
      list(
        id = "assay_tabs",
        width = 12,
        title = "Assays"
      ),
      lapply(assay_names, function(assay_name) {
        tabPanel(
          assay_name,
          assay_panel_ui(assay_name, obj)
        )
      })
    ))
  })
  
  observe({
    req(seurat_obj())
    obj <- seurat_obj()
    assay_names <- names(obj@assays)
    
    lapply(assay_names, function(assay_name) {
      assay_panel_server(assay_name, obj, config, output)
    })
  })
  
  output$reductions_ui <- renderUI({
    req(seurat_obj())
    obj <- seurat_obj()
    
    if (length(obj@reductions) == 0) {
      return(box(
        title = "Reductions",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        "No dimensional reductions found in this object."
      ))
    }
    
    reduction_names <- names(obj@reductions)
    
    fluidRow(
      column(12,
        box(
          title = "Dimensional Reductions",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          selectInput("selected_reduction", "Select Reduction:", 
                      choices = reduction_names,
                      selected = if(config$default_reduction %in% reduction_names) 
                        config$default_reduction else reduction_names[1]),
          selectInput("color_by", "Color by:", 
                      choices = c("None", colnames(obj@meta.data))),
          numericInput("dim1", "Dimension 1:", value = 1, min = 1),
          numericInput("dim2", "Dimension 2:", value = 2, min = 1),
          actionButton("plot_reduction", "Plot", class = "btn-primary")
        )
      ),
      column(12,
        box(
          title = "Reduction Plot",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          plotlyOutput("reduction_plot", height = "600px")
        )
      ),
      column(12,
        box(
          title = "Embeddings",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          DTOutput("embeddings_table")
        )
      )
    )
  })
  
  observeEvent(input$plot_reduction, {
    req(seurat_obj(), input$selected_reduction)
    
    output$reduction_plot <- renderPlotly({
      plot_reduction_interactive(
        seurat_obj(), 
        input$selected_reduction,
        input$color_by,
        c(input$dim1, input$dim2),
        config
      )
    })
    
    output$embeddings_table <- renderDT({
      obj <- seurat_obj()
      embeddings <- Embeddings(obj, reduction = input$selected_reduction)
      datatable(
        as.data.frame(embeddings),
        options = list(
          pageLength = config$rows_per_page,
          scrollX = TRUE
        )
      )
    })
  })
  
  output$metadata_ui <- renderUI({
    req(seurat_obj())
    obj <- seurat_obj()
    
    fluidRow(
      column(12,
        box(
          title = "Quick Navigation",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          div(style = "display: flex; gap: 15px; justify-content: center; padding: 10px;",
            actionButton("scroll_to_table", "Cell Metadata", 
                        icon = icon("table"),
                        class = "btn-primary",
                        style = "min-width: 150px;"),
            actionButton("scroll_to_summary", "Metadata Summary", 
                        icon = icon("chart-bar"),
                        class = "btn-info",
                        style = "min-width: 150px;"),
            actionButton("scroll_to_plots", "Distribution Plots", 
                        icon = icon("chart-area"),
                        class = "btn-success",
                        style = "min-width: 150px;")
          )
        )
      ),
      column(12,
        box(
          id = "metadata_table_box",
          title = "Cell Metadata",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          DTOutput("metadata_table")
        )
      ),
      column(12,
        box(
          id = "metadata_summary_box",
          title = "Metadata Summary",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          verbatimTextOutput("metadata_summary")
        )
      ),
      column(12,
        box(
          id = "metadata_plot_box",
          title = "Distribution Plots",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          selectInput("metadata_column", "Select Column:", 
                      choices = colnames(obj@meta.data)),
          plotOutput("metadata_plot", height = "400px")
        )
      )
    )
  })
  
  output$metadata_table <- renderDT({
    req(seurat_obj())
    obj <- seurat_obj()
    
    datatable(
      obj@meta.data,
      options = list(
        pageLength = config$rows_per_page,
        scrollX = TRUE,
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = '100px', targets = "_all"))
      ),
      filter = 'top',
      class = 'cell-border stripe'
    )
  })
  
  output$metadata_summary <- renderPrint({
    req(seurat_obj())
    summary(seurat_obj()@meta.data)
  })
  
  output$metadata_plot <- renderPlot({
    req(seurat_obj(), input$metadata_column)
    plot_metadata_distribution(seurat_obj(), input$metadata_column, config)
  })
  
  observeEvent(input$scroll_to_table, {
    shinyjs::runjs(
      "var el = document.getElementById('metadata_table_box'); 
       if(el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }"
    )
  })
  
  observeEvent(input$scroll_to_summary, {
    shinyjs::runjs(
      "var el = document.getElementById('metadata_summary_box'); 
       if(el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }"
    )
  })
  
  observeEvent(input$scroll_to_plots, {
    shinyjs::runjs(
      "var el = document.getElementById('metadata_plot_box'); 
       if(el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }"
    )
  })
  
  output$graphs_ui <- renderUI({
    req(seurat_obj())
    obj <- seurat_obj()
    
    if (length(obj@graphs) == 0) {
      return(box(
        title = "Graphs",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        "No graphs found in this object."
      ))
    }
    
    graph_names <- names(obj@graphs)
    
    fluidRow(
      column(12,
        box(
          title = "Neighbor Graphs",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          selectInput("selected_graph", "Select Graph:", choices = graph_names)
        )
      ),
      column(12,
        box(
          title = "Graph Information",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          verbatimTextOutput("graph_info")
        )
      )
    )
  })
  
  output$graph_info <- renderPrint({
    req(seurat_obj(), input$selected_graph)
    obj <- seurat_obj()
    graph <- obj@graphs[[input$selected_graph]]
    
    cat("Graph:", input$selected_graph, "\n")
    cat("Dimensions:", dim(graph), "\n")
    cat("Number of edges:", sum(graph > 0) / 2, "\n")
    cat("Sparsity:", 1 - (sum(graph > 0) / length(graph)), "\n")
  })
  
  output$images_ui <- renderUI({
    req(seurat_obj())
    obj <- seurat_obj()
    
    if (length(obj@images) == 0) {
      return(box(
        title = "Images",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        "No spatial images found in this object."
      ))
    }
    
    image_names <- names(obj@images)
    
    # Get all available features (genes) and metadata columns for coloring
    all_features <- rownames(obj)
    metadata_cols <- colnames(obj@meta.data)
    color_choices <- c("None", metadata_cols, all_features)
    
    fluidRow(
      column(12,
        box(
          title = "Spatial Images",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          selectInput("selected_image", "Select Image:", choices = image_names),
          selectInput("spatial_color_by", "Color by:", 
                      choices = color_choices,
                      selected = "None"),
          actionButton("plot_spatial", "Plot", class = "btn-primary"),
          hr(),
          plotOutput("spatial_plot", height = "600px")
        )
      )
    )
  })
  
  observeEvent(input$plot_spatial, {
    req(seurat_obj(), input$selected_image)
    
    output$spatial_plot <- renderPlot({
      obj <- seurat_obj()
      
      if (is.null(input$spatial_color_by) || input$spatial_color_by == "None") {
        SpatialPlot(obj, images = input$selected_image)
      } else {
        # Check if it's a metadata column or a feature
        if (input$spatial_color_by %in% colnames(obj@meta.data)) {
          SpatialPlot(obj, images = input$selected_image, group.by = input$spatial_color_by)
        } else {
          # It's a feature/gene
          SpatialFeaturePlot(obj, features = input$spatial_color_by, images = input$selected_image)
        }
      }
    })
  })
}

shinyApp(ui = ui, server = server)
