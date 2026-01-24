# Seurat Object Viewer

An interactive R Shiny dashboard for exploring Seurat v5 single-cell RNA-seq objects.

## Overview

This dashboard provides a comprehensive graphical interface to view all components of a Seurat object without performing additional analysis. It's designed for:
- Exploring processed scRNA-seq data
- Reviewing analysis results
- Sharing data with collaborators
- Quality control and validation

## Features

### 📊 Assay Views
- View `counts`, `data`, and `scale.data` matrices
- Feature and cell metadata
- Sparsity statistics
- Heatmaps of variable features

### 🎨 Reductions
- Interactive 2D plots (UMAP, tSNE, PCA, etc.)
- Color by any metadata column
- Overlay gene expression
- Export publication-ready plots

### 📋 Metadata
- Interactive cell metadata table
- Summary statistics
- Distribution plots
- Filter and search capabilities

### 🔗 Graphs
- Neighbor graph statistics
- Connectivity information
- Degree distributions

### 🗺️ Images (Spatial)
- Tissue image visualization
- Spatial feature expression
- Spot-level metadata overlay

## Installation

### Requirements
- R >= 4.0.0
- RStudio (recommended)

### R Packages
```r
install.packages(c("shiny", "shinydashboard", "DT", "ggplot2",
                   "plotly", "dplyr", "yaml", "qs2", "shinyjs", "Matrix"))

# Seurat v5
install.packages("Seurat")
```

## Usage

### Basic Usage
1. Start the app:
```r
# In R console
shiny::runApp("path/to/view_seurat")

# Or in RStudio, open app.R and click "Run App"

# Or use the ViewSeurat function
ViewSeurat()

# View a Seurat object directly
ViewSeurat(seurat_obj)

# View with custom title
ViewSeurat(seurat_obj, title = "My Analysis")
```

2. Upload your Seurat object (.rds or .qs2 file) via drag-and-drop or file browser
3. Explore the various components using the sidebar navigation

### Custom Configuration
1. Copy `inst/extdata/config.yaml.example` to `config.yaml` in the app directory
2. Edit `config.yaml` with your preferences:
   - Display limits
   - Color schemes
   - Plot parameters
   - Performance settings
3. Restart the app

### Example Configuration
```yaml
max_cells_display: 10000
default_reduction: "umap"
color_palette: "viridis"
point_size: 1.5
```

## File Requirements

### Input File
- **Formats**: 
  - `.rds` file containing a Seurat v5 object
  - `.qs2` file containing a Seurat v5 object (faster loading for large objects)
- **Creation**: 
  - RDS: `saveRDS(seurat_obj, "my_object.rds")`
  - QS2: `qs2::qs_save(seurat_obj, "my_object.qs2")`
- **Compatibility**: Must be a valid Seurat v5 object
- **Upload**: Drag and drop onto the upload zone or use the file browser

### Configuration File (Optional)
- **File**: `config.yaml` in the app directory
- **Template**: Copy from `inst/extdata/config.yaml.example`
- **Required**: No, app will use defaults if not present

## Troubleshooting

### Upload Issues
- Ensure file is a valid .rds or .qs2 file
- Check file size limits (default 10 GB)
- Verify Seurat object version compatibility
- For large objects (>1GB), consider using .qs2 format for faster loading

### Performance
- For large datasets (>50k cells), consider:
  - Reducing `max_cells_display` in config
  - Enabling sampling for matrix views
  - Using smaller subsets for initial exploration

### Visualization
- If plots are slow to render, reduce point size or enable downsampling
- For high-dimensional data, select specific features to display

## Project Structure
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

## Contributing
This is a visualization tool designed for exploring existing Seurat objects. If you need analysis capabilities, please use Seurat directly.

## License
MIT License

## Citation
If you use this tool, please cite:
- Seurat: Hao et al., Nature Biotechnology (2023)

## Support
For issues specific to:
- **Seurat objects**: See [Seurat documentation](https://satijalab.org/seurat/)
- **This viewer**: Check CLAUDE.md for technical details
