# viewseurat

Read-only Shiny app for viewing Seurat v4 and v5 objects. Upload an .rds or .qs2 file, or pass a Seurat object directly from R.

## Installation

```r
# install.packages("devtools")
devtools::install_github("pcantalupo/viewseurat")
```

## Usage

```r
# View an object directly
ViewSeurat(seurat_obj)

# Launch the upload interface (no object)
ViewSeurat()
```

You can also open `app.R` in RStudio and click "Run App" to get the upload interface.

## What you can view

- **Overview** — Object structure diagram, assay summary, and size on disk.
- **Assays** — counts, data, and scale.data matrices, variable features, feature metadata. Highlights the default assay.
- **Metadata** — Per-column summary with distribution graphics, searchable table, and distribution plots.
- **Reductions** — 2D scatter plots (UMAP, tSNE, PCA, etc.) colored by metadata. Shows which assay was used to build the reduction.
- **Images** — Spatial tissue plots colored by metadata. Supports 10X Visium (V1/V2) and FOV-based platforms (Nanostring CosMx, 10X Xenium).
- **The Guts** — Raw S4 slot inspection with clickable buttons for all slots

## Supported Seurat objects

- Seurat v4 and v5
- SCTAssay (sctransform)
- ChromatinAssay (Signac)
- Spatial- and image-based objects (10X Visium, Xenium, Nanostring CosMx)

This is a viewer only. It does not modify or analyze the Seurat object.
