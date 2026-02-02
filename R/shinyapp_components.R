#' ViewSeurat Shiny Application Components
#'
#' These functions create the UI and server components for the ViewSeurat Shiny app.
#'
#' @name viewseurat-app
#' @keywords internal
NULL

#' Create the ViewSeurat UI
#'
#' @return A Shiny UI object
#' @keywords internal
viewseurat_ui <- function() {

  shinydashboard::dashboardPage(
    skin = "blue",

    shinydashboard::dashboardHeader(
      title = "Seurat Object Viewer",
      shiny::tags$li(class = "dropdown",
              style = "padding: 15px; color: white; font-size: 16px;",
              shiny::uiOutput("file_title"))
    ),

    shinydashboard::dashboardSidebar(
      shinydashboard::sidebarMenu(id = "sidebar",
        shinydashboard::menuItem("Upload", tabName = "upload", icon = shiny::icon("upload")),
        shinydashboard::menuItem("Overview", tabName = "overview", icon = shiny::icon("info-circle")),
        shinydashboard::menuItem("Assays", tabName = "assays", icon = shiny::icon("table")),
        shinydashboard::menuItem("Reductions", tabName = "reductions", icon = shiny::icon("project-diagram")),
        shinydashboard::menuItem("Metadata", tabName = "metadata", icon = shiny::icon("list")),
        shinydashboard::menuItem("Graphs", tabName = "graphs", icon = shiny::icon("share-alt")),
        shinydashboard::menuItem("Images", tabName = "images", icon = shiny::icon("image")),
        shinydashboard::menuItem("The Guts", tabName = "guts", icon = shiny::icon("cogs"))
      )
    ),

    shinydashboard::dashboardBody(
      waiter::useWaiter(),
      shinyjs::useShinyjs(),
      shiny::tags$head(
        shiny::tags$style(shiny::HTML("
          .box-body { overflow-x: auto; }
          .dataTables_wrapper { overflow-x: auto; }
          .small-box { cursor: pointer; min-height: 120px; }
          .small-box h3 { font-size: 42px; }
          .small-box p { font-size: 16px; }
          .small-box .icon { font-size: 70px; }
          .shiny-file-input-progress {
            display: none;
          }
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
        shiny::tags$script(shiny::HTML("
          $(document).on('shiny:connected', function() {
            // Listen for file input change (both browse and drag-drop)
            $('#seurat_file').on('change', function(e) {
              if (this.files && this.files.length > 0) {
                // Immediately notify R to show waiter
                Shiny.setInputValue('file_upload_started', {
                  filename: this.files[0].name,
                  timestamp: Date.now()
                });
              }
            });

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
                  dropZone.innerHTML = '<h3><i class=\"fa fa-check-circle\"></i> File selected</h3><p>' + file.name + '</p>';
                } else {
                  alert('Please upload an .rds or .qs2 file');
                }
              }
            }
          });
        "))
      ),

      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "upload",
          shiny::fluidRow(
            shinydashboard::box(
              title = "Upload Seurat Object",
              status = "primary",
              solidHeader = TRUE,
              width = 12,
              shiny::div(
                id = "drop_zone",
                shiny::h3(shiny::icon("cloud-upload"), "Drag & Drop RDS or QS2 File Here"),
                shiny::p("or click below to browse"),
                shiny::br()
              ),
              shiny::fileInput("seurat_file",
                        "Choose RDS or QS2 file containing Seurat object",
                        accept = c(".rds", ".RDS", ".qs2", ".QS2")),
              shiny::verbatimTextOutput("upload_status")
            )
          )
        ),

        shinydashboard::tabItem(tabName = "overview",
          shiny::uiOutput("overview_ui")
        ),

        shinydashboard::tabItem(tabName = "assays",
          shiny::uiOutput("assays_ui")
        ),

        shinydashboard::tabItem(tabName = "reductions",
          shiny::uiOutput("reductions_ui")
        ),

        shinydashboard::tabItem(tabName = "metadata",
          shiny::uiOutput("metadata_ui")
        ),

        shinydashboard::tabItem(tabName = "graphs",
          shiny::uiOutput("graphs_ui")
        ),

        shinydashboard::tabItem(tabName = "images",
          shiny::uiOutput("images_ui")
        ),

        shinydashboard::tabItem(tabName = "guts",
          shiny::uiOutput("guts_ui")
        )
      )
    )
  )
}

#' Create the ViewSeurat Server Function
#'
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @keywords internal
viewseurat_server <- function(input, output, session) {

  seurat_obj <- shiny::reactiveVal(NULL)
  file_size_bytes <- shiny::reactiveVal(NULL)
  uploaded_filename <- shiny::reactiveVal("")

  # Check if object was pre-loaded via ViewSeurat()
  preloaded_obj <- shiny::getShinyOption("viewseurat.obj", default = NULL)
  preloaded_title <- shiny::getShinyOption("viewseurat.title", default = NULL)

  if (!is.null(preloaded_obj)) {
    seurat_obj(preloaded_obj)
    # Set title from parameter or default to "Pre-loaded Object"
    if (!is.null(preloaded_title) && preloaded_title != "") {
      uploaded_filename(preloaded_title)
    } else {
      uploaded_filename("")
    }
    # Auto-navigate to Overview tab when object is preloaded
    shinydashboard::updateTabItems(session, "sidebar", "overview")
  }

  output$file_title <- shiny::renderUI({
    if (uploaded_filename() != "") {
      shiny::tags$span(style = "color: white;", uploaded_filename())
    } else {
      NULL
    }
  })

  # Create waiter for upload progress
  upload_waiter <- waiter::Waiter$new(
    html = shiny::tagList(
      waiter::spin_fading_circles(),
      shiny::h4("Loading Seurat object..."),
      shiny::p("This may take several minutes for large files.")
    ),
    color = "rgba(0, 0, 0, 0.7)"
  )

  # Show waiter immediately when file selection detected
  shiny::observeEvent(input$file_upload_started, {
    upload_waiter$show()
  })

  # Process file when upload completes
  shiny::observeEvent(input$seurat_file, {
    shiny::req(input$seurat_file)

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

      # Store object
      seurat_obj(obj)
      file_size_bytes(input$seurat_file$size)
      uploaded_filename(input$seurat_file$name)

      output$upload_status <- shiny::renderText({
        paste0(
          "Successfully loaded Seurat object\n",
          "File type: ", toupper(file_ext), "\n",
          "Cells: ", ncol(obj), "\n",
          "Assays: ", paste(names(obj@assays), collapse = ", "), "\n",
          "Reductions: ", paste(names(obj@reductions), collapse = ", ")
        )
      })

      # Hide waiter and navigate
      upload_waiter$hide()
      shiny::showNotification("Seurat object loaded successfully!", type = "message")
      shinydashboard::updateTabItems(session, "sidebar", "overview")

    }, error = function(e) {
      upload_waiter$hide()
      output$upload_status <- shiny::renderText({
        paste("Error loading file:", e$message)
      })
      shiny::showNotification(paste("Error:", e$message), type = "error", duration = NULL)
    })
  })

  output$overview_ui <- shiny::renderUI({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    shiny::fluidRow(
      shiny::column(12,
        shinydashboard::box(
          title = "Standard Seurat summary",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          shiny::verbatimTextOutput("object_summary")
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          title = "View Seurat summary",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shiny::verbatimTextOutput("seurat_info_output")
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          title = "Size Info",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          shiny::verbatimTextOutput("size_info_output")
        )
      )
    )
  })

  output$object_summary <- shiny::renderPrint({
    shiny::req(seurat_obj())
    print(seurat_obj())
  })

  output$seurat_info_output <- shiny::renderPrint({
    shiny::req(seurat_obj())
    obj <- seurat_obj()
    info <- SeuratInfo(obj)

    cat(info$version, "\n")
    cat(info$graphs, "\n")
    cat(info$reductions, "\n")
    cat(info$images, "\n")
    cat(info$ident_label, "\n")
    cat("Idents():\n")
    print(info$idents_table)
    cat("\nAssays:\n")
    print(info$assays_table)
  })

  output$size_info_output <- shiny::renderPrint({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    # Get R session memory size
    mem_bytes <- as.numeric(object.size(obj))
    if (mem_bytes >= 1024^3) {
      mem_size <- paste0(round(mem_bytes / 1024^3, 2), " GB")
    } else if (mem_bytes >= 1024^2) {
      mem_size <- paste0(round(mem_bytes / 1024^2, 2), " MB")
    } else {
      mem_size <- paste0(round(mem_bytes / 1024, 2), " KB")
    }

    # Show file size on disk if available (from upload)
    disk_bytes <- file_size_bytes()
    if (!is.null(disk_bytes)) {
      if (disk_bytes >= 1024^3) {
        disk_size <- paste0(round(disk_bytes / 1024^3, 2), " GB")
      } else if (disk_bytes >= 1024^2) {
        disk_size <- paste0(round(disk_bytes / 1024^2, 2), " MB")
      } else {
        disk_size <- paste0(round(disk_bytes / 1024, 2), " KB")
      }
      cat("Size on disk:", disk_size, "\n")
    }

    cat("Size in R session:", mem_size, "\n")
  })

  shiny::observeEvent(input$goto_assays, {
    shinydashboard::updateTabItems(session, "sidebar", "assays")
  })

  shiny::observeEvent(input$goto_reductions, {
    shinydashboard::updateTabItems(session, "sidebar", "reductions")
  })

  shiny::observeEvent(input$goto_metadata, {
    shinydashboard::updateTabItems(session, "sidebar", "metadata")
  })

  output$assays_ui <- shiny::renderUI({
    shiny::req(seurat_obj())
    obj <- seurat_obj()
    assay_names <- names(obj@assays)
    default_assay <- Seurat::DefaultAssay(obj)

    # Create lightweight placeholder tabs - don't call assay_panel_ui() for all assays
    # Content is rendered lazily via uiOutput when tab is selected
    do.call(shinydashboard::tabBox, c(
      list(
        id = "assay_tabs",
        width = 12,
        title = "Assays"
      ),
      lapply(assay_names, function(assay_name) {
        # Add "(default)" indicator to the default assay tab (bolded)
        tab_label <- if (assay_name == default_assay) {
          shiny::tags$b(paste0(assay_name, " (default)"))
        } else {
          assay_name
        }
        shiny::tabPanel(
          tab_label,
          value = assay_name,  # Use assay name as value for programmatic selection
          # Use uiOutput placeholder instead of calling assay_panel_ui() directly
          shinycssloaders::withSpinner(
            shiny::uiOutput(paste0("assay_content_", assay_name))
          )
        )
      })
    ))
  })

  # Lazy initialization: only render assay panel UI and server when its tab is selected
  initialized_assays <- shiny::reactiveVal(character(0))

  # Helper function to initialize an assay (both UI and server)
  initialize_assay <- function(assay_name, obj) {
    if (!assay_name %in% initialized_assays()) {
      # Render the UI content for this assay
      output[[paste0("assay_content_", assay_name)]] <- shiny::renderUI({
        assay_panel_ui(assay_name, obj)
      })
      # Set up the server logic
      assay_panel_server(assay_name, obj, output)
      initialized_assays(c(initialized_assays(), assay_name))
    }
  }

  shiny::observeEvent(input$assay_tabs, {
    shiny::req(seurat_obj())
    obj <- seurat_obj()
    # Tab value is now the assay name directly
    selected_assay <- input$assay_tabs
    if (!is.null(selected_assay) && selected_assay %in% names(obj@assays)) {
      initialize_assay(selected_assay, obj)
    }
  })

  # Initialize and select the default assay when the Assays tab is first visited
  shiny::observeEvent(input$sidebar, {
    if (input$sidebar == "assays") {
      shiny::req(seurat_obj())
      obj <- seurat_obj()
      default_assay <- Seurat::DefaultAssay(obj)
      initialize_assay(default_assay, obj)
      # Auto-select the default assay tab
      shiny::updateTabsetPanel(session, "assay_tabs", selected = default_assay)
    }
  }, ignoreInit = TRUE)

  output$reductions_ui <- shiny::renderUI({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    if (length(obj@reductions) == 0) {
      return(shinydashboard::box(
        title = "Reductions",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        "No dimensional reductions found in this object."
      ))
    }

    reduction_names <- names(obj@reductions)

    # Get assay of origin for each reduction to display in dropdown
    # In selectInput, names are displayed and values are returned
    reduction_choices <- reduction_names
    names(reduction_choices) <- sapply(reduction_names, function(name) {
      assay_used <- obj[[name]]@assay.used
      paste0(name, " (", assay_used, ")")
    })

    shiny::fluidRow(
      shiny::column(12,
        shinydashboard::box(
          title = "Dimensional Reductions",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shiny::selectInput("selected_reduction", "Select Reduction:",
                      choices = reduction_choices,
                      selected = if("umap" %in% reduction_names) "umap" else reduction_names[1]),
          shiny::selectInput("color_by", "Color by:",
                      choices = c("None", colnames(obj@meta.data))),
          shiny::selectInput("dim1", "Dimension 1:", choices = 1, selected = 1),
          shiny::selectInput("dim2", "Dimension 2:", choices = 2, selected = 2),
          shiny::actionButton("plot_reduction", "Plot", class = "btn-primary")
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          title = "Reduction Plot",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          shinycssloaders::withSpinner(shiny::plotOutput("reduction_plot", height = "600px"))
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          title = "Embeddings",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          collapsed = TRUE,
          shinycssloaders::withSpinner(DT::DTOutput("embeddings_table"))
        )
      )
    )
  })

  # Update dimension dropdowns when reduction changes
  shiny::observeEvent(input$selected_reduction, {
    shiny::req(seurat_obj(), input$selected_reduction)
    obj <- seurat_obj()
    embeddings <- Seurat::Embeddings(obj, reduction = input$selected_reduction)
    n_dims <- ncol(embeddings)
    dim_choices <- as.character(1:n_dims)

    shiny::updateSelectInput(session, "dim1", choices = dim_choices, selected = "1")
    shiny::updateSelectInput(session, "dim2", choices = dim_choices, selected = if(n_dims >= 2) "2" else "1")
  })

  shiny::observeEvent(input$plot_reduction, {
    shiny::req(seurat_obj(), input$selected_reduction)

    output$reduction_plot <- shiny::renderPlot({
      plot_reduction_static(
        seurat_obj(),
        input$selected_reduction,
        input$color_by,
        c(as.integer(input$dim1), as.integer(input$dim2))
      )
    })

    output$embeddings_table <- DT::renderDT({
      obj <- seurat_obj()
      embeddings <- Seurat::Embeddings(obj, reduction = input$selected_reduction)
      DT::datatable(
        as.data.frame(embeddings),
        options = list(
          pageLength = 10,
          scrollX = TRUE
        )
      )
    }, server = TRUE)
  })

  output$metadata_ui <- shiny::renderUI({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    shiny::fluidRow(
      shiny::column(12,
        shinydashboard::box(
          title = "Quick Navigation",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shiny::div(style = "display: flex; gap: 15px; justify-content: center; padding: 10px;",
            shiny::actionButton("scroll_to_table", "Cell Metadata",
                        icon = shiny::icon("table"),
                        class = "btn-primary",
                        style = "min-width: 150px;"),
            shiny::actionButton("scroll_to_summary", "Metadata Summary",
                        icon = shiny::icon("chart-bar"),
                        class = "btn-info",
                        style = "min-width: 150px;"),
            shiny::actionButton("scroll_to_plots", "Distribution Plots",
                        icon = shiny::icon("chart-area"),
                        class = "btn-success",
                        style = "min-width: 150px;")
          )
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          id = "metadata_table_box",
          title = "Cell Metadata",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shinycssloaders::withSpinner(DT::DTOutput("metadata_table"))
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          id = "metadata_summary_box",
          title = "Metadata Summary",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          shiny::verbatimTextOutput("metadata_summary")
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          id = "metadata_plot_box",
          title = "Distribution Plots",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          shiny::selectInput("metadata_column", "Select Column:",
                      choices = colnames(obj@meta.data)),
          shiny::plotOutput("metadata_plot", height = "400px")
        )
      )
    )
  })

  output$metadata_table <- DT::renderDT({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    DT::datatable(
      obj@meta.data,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = '100px', targets = "_all"))
      ),
      filter = 'top',
      class = 'cell-border stripe'
    )
  }, server = TRUE)

  output$metadata_summary <- shiny::renderPrint({
    shiny::req(seurat_obj())
    summary(seurat_obj()@meta.data)
  })

  output$metadata_plot <- shiny::renderPlot({
    shiny::req(seurat_obj(), input$metadata_column)
    plot_metadata_distribution(seurat_obj(), input$metadata_column)
  })

  shiny::observeEvent(input$scroll_to_table, {
    shinyjs::runjs(
      "var el = document.getElementById('metadata_table_box');
       if(el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }"
    )
  })

  shiny::observeEvent(input$scroll_to_summary, {
    shinyjs::runjs(
      "var el = document.getElementById('metadata_summary_box');
       if(el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }"
    )
  })

  shiny::observeEvent(input$scroll_to_plots, {
    shinyjs::runjs(
      "var el = document.getElementById('metadata_plot_box');
       if(el) { el.scrollIntoView({behavior: 'smooth', block: 'start'}); }"
    )
  })

  output$graphs_ui <- shiny::renderUI({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    if (length(obj@graphs) == 0) {
      return(shinydashboard::box(
        title = "Graphs",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        "No graphs found in this object."
      ))
    }

    graph_names <- names(obj@graphs)

    shiny::fluidRow(
      shiny::column(12,
        shinydashboard::box(
          title = "Neighbor Graphs",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shiny::selectInput("selected_graph", "Select Graph:", choices = graph_names)
        )
      ),
      shiny::column(12,
        shinydashboard::box(
          title = "Graph Information",
          status = "info",
          solidHeader = TRUE,
          width = 12,
          shiny::verbatimTextOutput("graph_info")
        )
      )
    )
  })

  output$graph_info <- shiny::renderPrint({
    shiny::req(seurat_obj(), input$selected_graph)
    obj <- seurat_obj()
    graph <- obj@graphs[[input$selected_graph]]

    cat("Graph:", input$selected_graph, "\n")
    cat("Dimensions:", dim(graph), "\n")
    cat("Number of edges:", sum(graph > 0) / 2, "\n")
    cat("Sparsity:", 1 - (sum(graph > 0) / length(graph)), "\n")
  })

  output$images_ui <- shiny::renderUI({
    shiny::req(seurat_obj())
    obj <- seurat_obj()

    if (length(obj@images) == 0) {
      return(shinydashboard::box(
        title = "Images",
        status = "warning",
        solidHeader = TRUE,
        width = 12,
        "No spatial images found in this object."
      ))
    }

    image_names <- names(obj@images)

    # Only get metadata columns initially - genes loaded via server-side selectize
    metadata_cols <- colnames(obj@meta.data)

    shiny::fluidRow(
      shiny::column(12,
        shinydashboard::box(
          title = "Spatial Images",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shiny::selectInput("selected_image", "Select Image:", choices = image_names),
          shiny::selectizeInput("spatial_color_by", "Color by:",
                      choices = c("None", metadata_cols),
                      selected = "None",
                      options = list(
                        placeholder = "Select or type gene name..."
                      )),
          shiny::actionButton("plot_spatial", "Plot", class = "btn-primary"),
          shiny::hr(),
          shiny::plotOutput("spatial_plot", height = "600px")
        )
      )
    )
  })

  # Load gene choices in background via server-side selectize
  shiny::observeEvent(seurat_obj(), {
    obj <- seurat_obj()
    if (length(obj@images) > 0) {
      all_choices <- c("None", colnames(obj@meta.data), rownames(obj))
      shiny::updateSelectizeInput(session, "spatial_color_by",
        choices = all_choices,
        selected = "None",
        server = TRUE)
    }
  })

  shiny::observeEvent(input$plot_spatial, {
    shiny::req(seurat_obj(), input$selected_image)

    output$spatial_plot <- shiny::renderPlot({
      obj <- seurat_obj()

      if (is.null(input$spatial_color_by) || input$spatial_color_by == "None") {
        Seurat::SpatialPlot(obj, images = input$selected_image)
      } else if (input$spatial_color_by %in% colnames(obj@meta.data)) {
        # Metadata column - check if numeric or categorical
        if (is.numeric(obj@meta.data[[input$spatial_color_by]])) {
          # Numeric metadata: use SpatialFeaturePlot for continuous color scale
          Seurat::SpatialFeaturePlot(obj, features = input$spatial_color_by, images = input$selected_image)
        } else {
          # Categorical metadata: use SpatialPlot with group.by
          Seurat::SpatialPlot(obj, images = input$selected_image, group.by = input$spatial_color_by)
        }
      } else {
        # It's a feature/gene
        Seurat::SpatialFeaturePlot(obj, features = input$spatial_color_by, images = input$selected_image)
      }
    })
  })

  output$guts_ui <- shiny::renderUI({
    shiny::req(seurat_obj())

    shiny::fluidRow(
      shiny::column(12,
        shinydashboard::box(
          title = "Object Structure",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          shiny::p("Internal structure of the Seurat object using str():"),
          shiny::verbatimTextOutput("guts_output")
        )
      )
    )
  })

  output$guts_output <- shiny::renderPrint({
    shiny::req(seurat_obj())
    str(seurat_obj())
  })
}

#' Create the ViewSeurat Shiny Application
#'
#' This function creates and returns the Shiny app object for viewing Seurat objects.
#' It is called internally by \code{\link{ViewSeurat}}.
#'
#' @return A Shiny app object
#' @keywords internal
viewseurat_app <- function() {
  # Set max upload size (10 GB)
  options(shiny.maxRequestSize = 10240 * 1024^2)

  shiny::shinyApp(ui = viewseurat_ui(), server = viewseurat_server)
}
