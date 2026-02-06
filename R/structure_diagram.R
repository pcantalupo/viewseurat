#' Structure Diagram UI
#'
#' Generate a color-coded HTML/CSS block diagram of a Seurat object's structure.
#' Shows assays with their layers and dimensions, metadata, reductions, graphs,
#' and spatial images in a visual layout similar to Bioconductor schematics.
#'
#' @param obj A Seurat object
#' @return A \code{shiny::tagList} containing the rendered diagram
#' @keywords internal
structure_diagram_ui <- function(obj) {
  n_cells <- ncol(obj)
  n_features <- nrow(obj)
  assay_names <- names(obj@assays)
  default_assay <- Seurat::DefaultAssay(obj)
  reduction_names <- names(obj@reductions)
  graph_names <- names(obj@graphs)
  image_names <- names(obj@images)
  n_meta_cols <- ncol(obj@meta.data)

  # --- CSS ---
  css <- shiny::tags$style(shiny::HTML("
    .vs-structure-outer {
      border: 2px solid #bdbdbd;
      border-radius: 10px;
      background: #fafafa;
      padding: 20px;
      max-width: 900px;
      font-family: 'Helvetica Neue', Arial, sans-serif;
    }
    .vs-structure-title {
      font-size: 20px;
      font-weight: 700;
      margin-bottom: 2px;
    }
    .vs-structure-subtitle {
      font-size: 14px;
      color: #616161;
      margin-bottom: 18px;
    }
    .vs-section {
      border-radius: 8px;
      padding: 12px 14px;
      margin-bottom: 14px;
    }
    .vs-section-label {
      font-weight: 700;
      font-size: 13px;
      margin-bottom: 8px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }
    .vs-assays-section {
      background: #e3f2fd;
      border: 1.5px solid #1976d2;
    }
    .vs-assays-section .vs-section-label { color: #1976d2; }
    .vs-assay-cards {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
    }
    .vs-assay-card {
      background: #fff;
      border-left: 4px solid #1976d2;
      border-radius: 6px;
      padding: 10px 14px;
      min-width: 140px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.08);
    }
    .vs-assay-card.vs-default {
      border-left-color: #0d47a1;
      background: #e8eaf6;
    }
    .vs-assay-name {
      font-weight: 700;
      font-size: 14px;
      color: #1565c0;
    }
    .vs-assay-card.vs-default .vs-assay-name { color: #0d47a1; }
    .vs-assay-dim {
      font-size: 12px;
      color: #616161;
      margin: 2px 0 4px 0;
    }
    .vs-assay-layers {
      font-size: 11px;
      color: #757575;
    }
    .vs-assay-layer-item {
      display: inline-block;
      background: #e3f2fd;
      border-radius: 3px;
      padding: 1px 6px;
      margin: 2px 2px 0 0;
      font-size: 11px;
      color: #1565c0;
    }
    .vs-bottom-row {
      display: flex;
      flex-wrap: wrap;
      gap: 12px;
    }
    .vs-bottom-row > .vs-section {
      flex: 1 1 200px;
      min-width: 180px;
      margin-bottom: 0;
    }
    .vs-meta-section {
      background: #e8f5e9;
      border: 1.5px solid #388e3c;
    }
    .vs-meta-section .vs-section-label { color: #388e3c; }
    .vs-reductions-section {
      background: #fff3e0;
      border: 1.5px solid #f57c00;
    }
    .vs-reductions-section .vs-section-label { color: #f57c00; }
    .vs-graphs-section {
      background: #f3e5f5;
      border: 1.5px solid #7b1fa2;
    }
    .vs-graphs-section .vs-section-label { color: #7b1fa2; }
    .vs-images-section {
      background: #fce4ec;
      border: 1.5px solid #c2185b;
    }
    .vs-images-section .vs-section-label { color: #c2185b; }
    .vs-item {
      font-size: 13px;
      color: #424242;
    }
    .vs-empty {
      font-size: 13px;
      color: #9e9e9e;
      font-style: italic;
    }
    .vs-badge {
      display: inline-block;
      background: rgba(0,0,0,0.07);
      border-radius: 3px;
      padding: 1px 6px;
      margin: 2px 3px 0 0;
      font-size: 12px;
    }
  "))

  # --- Assay cards ---
  assay_cards <- lapply(assay_names, function(aname) {
    assay_obj <- obj@assays[[aname]]
    layers <- SeuratObject::Layers(assay_obj)
    is_default <- aname == default_assay

    # Get assay-level dimensions
    a_nrow <- nrow(assay_obj)
    a_ncol <- ncol(assay_obj)

    # Per-layer dimension tags
    layer_tags <- lapply(layers, function(lname) {
      ld <- get_layer_dim(assay_obj, lname)
      label <- if (ld[1] > 0 && ld[2] > 0) {
        paste0(lname, " [", format_number(ld[1]), " x ", format_number(ld[2]), "]")
      } else {
        lname
      }
      shiny::tags$span(class = "vs-assay-layer-item", label)
    })

    card_class <- if (is_default) "vs-assay-card vs-default" else "vs-assay-card"
    name_label <- if (is_default) paste0(aname, " *") else aname

    shiny::tags$div(class = card_class,
      shiny::tags$div(class = "vs-assay-name", name_label),
      shiny::tags$div(class = "vs-assay-dim",
        paste0(format_number(a_nrow), " features x ", format_number(a_ncol), " cells")
      ),
      shiny::tags$div(class = "vs-assay-layers", layer_tags)
    )
  })

  assays_section <- shiny::tags$div(class = "vs-section vs-assays-section",
    shiny::tags$div(class = "vs-section-label", "Assays"),
    shiny::tags$div(class = "vs-assay-cards", assay_cards)
  )

  # --- Metadata section ---
  meta_section <- shiny::tags$div(class = "vs-section vs-meta-section",
    shiny::tags$div(class = "vs-section-label", "Metadata"),
    shiny::tags$div(class = "vs-item",
      paste0(format_number(n_cells), " cells"),
      shiny::tags$br(),
      paste0(n_meta_cols, " columns")
    )
  )

  # --- Reductions section ---
  if (length(reduction_names) > 0) {
    red_badges <- lapply(reduction_names, function(rname) {
      n_dims <- ncol(obj@reductions[[rname]])
      shiny::tags$span(class = "vs-badge", paste0(rname, " (", n_dims, ")"))
    })
    red_content <- shiny::tags$div(class = "vs-item", red_badges)
  } else {
    red_content <- shiny::tags$div(class = "vs-empty", "None")
  }
  reductions_section <- shiny::tags$div(class = "vs-section vs-reductions-section",
    shiny::tags$div(class = "vs-section-label", "Reductions"),
    red_content
  )

  # --- Graphs section ---
  if (length(graph_names) > 0) {
    graph_badges <- lapply(graph_names, function(gname) {
      shiny::tags$span(class = "vs-badge", gname)
    })
    graph_content <- shiny::tags$div(class = "vs-item", graph_badges)
  } else {
    graph_content <- shiny::tags$div(class = "vs-empty", "None")
  }
  graphs_section <- shiny::tags$div(class = "vs-section vs-graphs-section",
    shiny::tags$div(class = "vs-section-label", "Graphs"),
    graph_content
  )

  # --- Images section ---
  if (length(image_names) > 0) {
    img_badges <- lapply(image_names, function(iname) {
      shiny::tags$span(class = "vs-badge", iname)
    })
    img_content <- shiny::tags$div(class = "vs-item", img_badges)
  } else {
    img_content <- shiny::tags$div(class = "vs-empty", "None")
  }
  images_section <- shiny::tags$div(class = "vs-section vs-images-section",
    shiny::tags$div(class = "vs-section-label", "Images"),
    img_content
  )

  # --- Assemble ---
  shiny::tagList(
    css,
    shiny::tags$div(class = "vs-structure-outer",
      shiny::tags$div(class = "vs-structure-title", "Seurat Object"),
      shiny::tags$div(class = "vs-structure-subtitle",
        paste0(
          format_number(n_features), " features x ",
          format_number(n_cells), " cells",
          "    |    * = default assay"
        )
      ),
      assays_section,
      shiny::tags$div(class = "vs-bottom-row",
        meta_section,
        reductions_section
      ),
      shiny::tags$div(class = "vs-bottom-row",
        graphs_section,
        images_section
      )
    )
  )
}
