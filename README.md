# Seurat Object Viewer

An interactive R Shiny dashboard for exploring Seurat v5 single-cell RNA-seq objects.

# Overview

This dashboard provides a comprehensive graphical interface to view all components of a Seurat object without performing additional analysis. It's designed for:
- Exploring the internal structure of Seurat objects
- Reviewing metadata and analysis results
- Quality control and validation

# Features

## 🧬 Overview
- Summary of object structure and contents
- Object-level metadata and statistics

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

# Miscellaneous
This is a visualization tool designed for exploring existing Seurat objects. If you need analysis capabilities, please use the Seurat package directly.


