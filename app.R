# Shiny App Launcher
#
# This file supports two development workflows:
#
# 1. RStudio "Run App" button: Click "Run App" to launch the upload interface
#    Note: A warning about "Loading R/ subdirectory" is expected and harmless -
#    we intentionally source R/ files manually below.
# 2. Package development: Use devtools::load_all() then ViewSeurat(obj)
#
# For production use after installing the package:
#   viewseurat::ViewSeurat()
#   viewseurat::ViewSeurat(seurat_obj)
#   viewseurat::ViewSeurat(seurat_obj, title = "My Analysis")

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

# Source all package R files
source("R/seurat_utils.R")
source("R/plot_functions.R")
source("R/assay_panel.R")
source("R/shinyapp_components.R")

# Set max upload size (10 GB)
options(shiny.maxRequestSize = 10240 * 1024^2)

# Create and return the app (explicit shinyApp call for RStudio detection)
shinyApp(ui = viewseurat_ui(), server = viewseurat_server)
