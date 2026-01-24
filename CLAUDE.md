# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Shiny dashboard application for viewing Seurat v5 single-cell RNA-seq objects. It is a read-only viewer with NO analysis capabilities - it only displays existing data from processed Seurat objects.

## Running the Application

There are two ways to run the app:

### Method 1: RStudio "Run App" (Easy usage without installation)
```r
# Open app.R in RStudio and click "Run App"
# This sources all R/ files directly and launches shinyApp()
```

### Method 2: ViewSeurat() Function (Interactive use / Development)
```r
# Install the package first
devtools::install_local(".")

# Then use the function
ViewSeurat(seurat_obj)  # Direct viewing
ViewSeurat(seurat_obj, title = "My Analysis")  # Custom title
ViewSeurat()  # Launch app with file upload prompt
```

The app accepts `.rds` or `.qs2` files containing Seurat v5 objects via drag-and-drop upload.

### File Organization (Mastering Shiny Pattern)

This package follows the pattern from Hadley Wickham's "Mastering Shiny" book:

- **`app.R` (root)**: Manual launcher - sources R/ files and calls `shinyApp(ui, server)`
- **`R/shinyapp_components.R`**: Contains `viewseurat_ui()` and `viewseurat_server()` that define the Shiny app
- **`R/ViewSeurat.R`**: Contains `ViewSeurat()` exported function that wraps the app

When making changes to the app, **edit `R/shinyapp_components.R`** since it contains the ui and server logic.

## Architecture

### Application Structure

The app follows a modular Shiny architecture with the app wrapped in a function:
- `R/shinyapp_components.R` - `viewseurat_ui()` and `viewseurat_server()` defining the Shiny app
- `R/seurat_utils.R` - Utilities for extracting and processing Seurat object data
- `R/plot_functions.R` - Plotting functions for visualizations
- `R/ui_modules.R` - Reusable UI module functions for assay panels
- `R/ViewSeurat.R` - `ViewSeurat()` exported function that wraps the app

### Configuration System

Configuration is managed through YAML:
- `inst/extdata/config.yaml.example` - Template with all available options and defaults
- `load_config()` in seurat_utils.R loads config with fallback to defaults
- Config controls: display limits, plot aesthetics, performance settings, file upload size

Access config in server functions via the `config` variable passed as a parameter.

The configuration system is not fully implemented for user overrides yet.

### Key Reactive Flow

1. User uploads file → `input$seurat_file` observeEvent triggers
2. File loaded and validated → stored in `seurat_obj` reactiveVal
3. UI sections render dynamically based on seurat_obj() content
4. Each tab (Overview, Assays, Reductions, Metadata, Graphs, Images) uses `renderUI()` to generate UI on demand
5. Assay panels use module pattern: `assay_panel_ui()` and `assay_panel_server()` for each assay

### Seurat Object Structure

The viewer displays these Seurat v5 components:
- **Assays** (`obj@assays`): Expression matrices with layers (counts, data, scale.data), variable features, and feature metadata
- **Reductions** (`obj@reductions`): Dimensional reductions (UMAP, tSNE, PCA) with embeddings
- **Metadata** (`obj@meta.data`): Cell-level annotations and QC metrics
- **Graphs** (`obj@graphs`): Neighbor graphs for clustering
- **Images** (`obj@images`): Spatial transcriptomics tissue images


### Working with Spatial Data

Spatial plots use Seurat's built-in functions:
- `SpatialPlot()` for basic tissue visualization with metadata overlay
- `SpatialFeaturePlot()` for gene expression overlay on tissue images
- Check `input$spatial_color_by` against metadata columns vs. feature names to determine which function to use

### Data Layer Access (Seurat v5)

Seurat v5 uses a `layers` slot instead of direct slots. Access data using:
```r
GetAssayData(obj, layer = "counts")     # Raw counts, or "data", "scale.data"
```

The `get_assay_data_safe()` function wraps this with error handling.

### Sparse Matrix Handling

Most Seurat data uses sparse matrices (dgCMatrix):
- Use `get_matrix_sample()` to extract small representative samples for display
- Use `get_sparsity_info()` to calculate memory and sparsity statistics
- Convert to dense matrix only for small samples: `as.matrix(sparse_matrix[rows, cols])`

## Common Development Patterns

### Plot Rendering Strategy

Two rendering approaches:
1. **Static plots** (ggplot2): Used for metadata distributions and heatmaps
2. **Interactive plots** (plotly): Used for dimensional reductions via `plot_reduction_interactive()`

For large datasets, sampling is applied automatically based on config limits.

### Adding a New Viewer Tab

1. Add menuItem to sidebar in `dashboardSidebar()` in `R/shinyapp_components.R`
2. Add tabItem to `tabItems()` in dashboardBody in `R/shinyapp_components.R`
3. Create corresponding `output$[name]_ui <- renderUI({})` in the server function
4. Implement UI generation logic that checks for data availability

### Handling Missing Data

Always check for component existence before rendering:
```r
if (length(obj@reductions) == 0) {
  return(shinydashboard::box(title = "Reductions", status = "warning",
             "No dimensional reductions found"))
}
```

## Miscellaneous Details
### File Upload Handling

- Drag-and-drop implemented via JavaScript in `R/shinyapp_components.R`
- Accepts `.rds` (via `readRDS()`) and `.qs2` (via `qs2::qs_read()`)
- Max file size controlled by `options(shiny.maxRequestSize = config$max_upload_size_mb * 1024^2)`
- Validation: Check file extension, verify object is class "Seurat"

### Performance Considerations

- Default max upload: 10GB (configurable)
- Matrix display: Sample 50 rows x 20 cols by default (configurable)
- Reduction plots: Use plotly for efficient rendering of large point clouds
- All DataTables use server-side pagination with configurable `rows_per_page`

### Important Constraints

- **Read-only viewer**: Do NOT add analysis functions, data transformations, or object modifications
- **Seurat v5 only**: Code assumes v5 object structure (layers vs. direct slots)
- **No external data**: App only displays what's in the uploaded Seurat object
- **No persistence**: Uploaded objects are session-only (not saved server-side)
