# Seurat Object Viewer

An interactive R Shiny dashboard for exploring Seurat v5 single-cell RNA-seq objects.

# Overview

This dashboard provides a comprehensive graphical interface to view all components of a Seurat object without performing additional analysis. It's designed for:
- Exploring the internal structure of Seurat objects
- Reviewing analysis results
- Quality control and validation

# Features

## 📊 Assays
- At a glance summary of all assays
- View the `counts`, `data`, and `scale.data` matrices
- Feature and cell metadata
- Sparsity statistics

## 🎨 Reductions
- Interactive 2D plots (UMAP, tSNE, PCA, etc.)
- Color by metadata columns

## 📋 Metadata
- Interactive cell metadata table with search capabilities
- Summary statistics and distribution plots

## 🔗 Graphs
- Neighbor graph statistics (dimensions, edges, sparsity)

## 🗺️ Images (Spatial)
- Tissue image visualization
- Spot-level metadata overlay

# Installation
Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("pcantalupo/viewseurat")
```

# Usage

## Basic Usage
1. Start the app:
```r
# In R console
shiny::runApp("path/to/viewseurat")

# Or in RStudio, open app.R and click "Run App"

# Or use the ViewSeurat() to view a Seurat object directly
ViewSeurat(seurat_obj)

# View with custom title
ViewSeurat(seurat_obj, title = "My Seurat Object")
```

2. Upload your Seurat object (.rds or .qs2 file) via drag-and-drop or file browser
3. Explore the various components using the sidebar navigation

## Custom Configuration
**Under Development**: Configuration options will be expanded in future releases.

1. Copy `inst/extdata/config.yaml.example` to `config.yaml` in the app directory
2. Edit `config.yaml` with your preferences:
   - Display limits
   - Color schemes
   - Plot parameters
   - Performance settings
3. Restart the app

# Project Structure
```
viewseurat/
├── app.R                           # Development launcher (sources R/, runs shinyApp)
├── DESCRIPTION                     # Package metadata
├── R/
│   ├── shinyapp_components.R       # Main UI and server logic
│   ├── ViewSeurat.R                # ViewSeurat() exported function
│   ├── seurat_utils.R              # Seurat object utilities
│   ├── plot_functions.R            # Plotting functions
│   └── ui_modules.R                # Shiny UI modules
└── inst/
    └── extdata/
        └── config.yaml.example     # Example configuration
```

# Miscellaneous
This is a visualization tool designed for exploring existing Seurat objects. If you need analysis capabilities, please use the Seurat package directly.


