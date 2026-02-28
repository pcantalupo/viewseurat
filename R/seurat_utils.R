#' Get Top Variable Features
#'
#' Get the top variable features for an assay.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay
#' @param n Number of features to return
#' @return Character vector of feature names
#' @keywords internal
#' @export
#' @importFrom Seurat VariableFeatures DefaultAssay<-
get_top_variable_features <- function(obj, assay_name, n = 20) {
  DefaultAssay(obj) <- assay_name
  var_features <- VariableFeatures(obj)

  if (length(var_features) == 0) {
    return(NULL)
  }

  head(var_features, n)
}

#' Format Number
#'
#' Format large numbers with comma-separated thousands.
#'
#' @param x A numeric value
#' @return A formatted string
#' @keywords internal
#' @export
format_number <- function(x) {
  format(as.integer(x), big.mark = ",")
}

#' Validate Seurat Object
#'
#' Check that a Seurat object is valid for viewing.
#'
#' @param obj An object to validate
#' @return TRUE if valid, otherwise throws an error
#' @keywords internal
#' @export
validate_seurat_object <- function(obj) {
  # Check type first and return early - subsequent checks would fail with cryptic errors
 if (!inherits(obj, "Seurat")) {
    stop("Object is not a Seurat object")
  }

  errors <- c()

  # Must check assays first: ncol() dispatches through [[DefaultAssay()]] which
  # throws a cryptic error when @assays is empty.
  if (length(obj@assays) == 0) {
    errors <- c(errors, "Object contains no assays")
  } else if (ncol(obj) == 0) {
    errors <- c(errors, "Object contains no cells")
  }

  if (length(errors) > 0) {
    stop(paste(errors, collapse = "; "))
  }

  return(TRUE)
}

#' Get Assay Data Safely with Subsetting
#'
#' Wrapper around LayerData with error handling and optional subsetting.
#' Subsetting at retrieval time is much faster for large/on-disk matrices
#' because LayerData() with features/cells params avoids loading the full matrix.
#'
#' @param obj A Seurat object
#' @param assay_name Name of the assay
#' @param layer Name of the layer (counts, data, scale.data)
#' @param max_features Maximum features to retrieve (NULL for all)
#' @param max_cells Maximum cells to retrieve (NULL for all)
#' @param cells_filter Integer vector of cell indices to retrieve; overrides
#'   max_cells when provided.
#' @return Matrix data or NULL if not available
#' @keywords internal
#' @export
#' @importFrom SeuratObject LayerData
get_assay_data_safe <- function(obj, assay_name, layer,
                                max_features = NULL, max_cells = NULL,
                                cells_filter = NULL) {
  tryCatch({
    assay <- obj@assays[[assay_name]]
    if (is.null(assay)) return(NULL)

    # Determine subsetting indices
    features <- NULL
    cells <- NULL

    if (!is.null(max_features)) {
      n_features <- nrow(assay)
      if (max_features < n_features) {
        features <- 1:max_features
      }
    }

    if (!is.null(cells_filter)) {
      cells <- cells_filter
    } else if (!is.null(max_cells)) {
      n_cells <- ncol(assay)
      if (max_cells < n_cells) {
        cells <- 1:max_cells
      }
    }

    # Old-style Assay subclasses (e.g. Signac ChromatinAssay) lack a @layers
    # slot and their LayerData S3 method may silently ignore features/cells
    # params, returning the full matrix -> OOM. Bypass LayerData entirely and
    # use direct slot access + manual subsetting for these assay types.
    use_slot_fallback <- !"layers" %in% methods::slotNames(assay) &&
      !methods::is(assay, "Assay5")

    if (use_slot_fallback) {
      if (!layer %in% methods::slotNames(assay)) return(NULL)
      mat <- methods::slot(assay, layer)
      if (is.null(mat) || length(dim(mat)) < 2) return(NULL)
      row_idx <- if (!is.null(features)) features else seq_len(nrow(mat))
      col_idx <- if (!is.null(cells)) cells else seq_len(ncol(mat))
      mat[row_idx, col_idx, drop = FALSE]
    } else {
      # LayerData returns an empty matrix (not NULL/error) for missing layers.
      # Guard explicitly so callers get a clean NULL.
      if (!layer %in% SeuratObject::Layers(assay)) return(NULL)
      # LayerData subsets at retrieval - much faster for large/on-disk data
      SeuratObject::LayerData(assay, layer = layer,
                              features = features, cells = cells)
    }
  }, error = function(e) {
    NULL
  })
}

#' Find Ident Label
#'
#' Find the metadata column name that matches the active Idents.
#'
#' @param seurat A Seurat object
#' @return The name of the matching metadata column, or "Unknown"
#' @keywords internal
#' @export
#' @importFrom Seurat Idents
FindIdentLabel <- function(seurat) {
  # Find the metadata column name that matches Idents
  ident_values <- as.character(Idents(seurat))

  for (col in colnames(seurat@meta.data)) {
    if (identical(as.character(seurat@meta.data[[col]]), ident_values)) {
      return(col)
    }
  }

  return("Unknown")
}

#' Get Layer Dimensions
#'
#' Returns the dimensions of a specific layer within an assay, handling both
#' v5 Assay objects (layers slot) and v4-style direct slots. Works on SCTAssay
#' and other assay subclasses that may lack a layers slot.
#'
#' @param assay_obj A Seurat Assay object
#' @param layer Name of the layer (e.g. "counts", "data", "scale.data")
#' @return Integer vector c(nrow, ncol), or c(0L, 0L) if layer is absent
#' @keywords internal
#' @export
get_layer_dim <- function(assay_obj, layer) {
  slotnames <- methods::slotNames(assay_obj)
  mat <- NULL
  if ("layers" %in% slotnames) {
    mat <- assay_obj@layers[[layer]]
  } else if (layer %in% slotnames) {
    mat <- methods::slot(assay_obj, layer)
  }

  if (is.null(mat)) return(c(0L, 0L))
  dim(mat)
}

#' Build Metadata Column Profile
#'
#' Compute a summary data frame with one row per metadata column, including
#' type, inline HTML distribution bar, completeness, and descriptive stats.
#'
#' @param meta A data.frame (typically obj@@meta.data)
#' @param max_levels Maximum number of factor/character levels to display
#' @return A data.frame suitable for DT::datatable with escape = FALSE
#' @keywords internal
#' @export
build_metadata_profile <- function(meta, max_levels = 5L) {
  n_rows <- nrow(meta)

  profiles <- lapply(colnames(meta), function(col) {
    vals <- meta[[col]]
    col_class <- class(vals)[1]
    n_na <- sum(is.na(vals))
    pct_complete <- round((1 - n_na / n_rows) * 100, 1)

    if (is.numeric(vals)) {
      info <- .profile_numeric(vals, n_na, n_rows)
      type_label <- col_class
    } else if (is.logical(vals)) {
      info <- .profile_logical(vals, n_na, n_rows)
      type_label <- "logical"
    } else if (is.factor(vals)) {
      info <- .profile_categorical(vals, n_rows, max_levels)
      type_label <- "factor"
    } else {
      # character or other -- treat as factor
      info <- .profile_categorical(vals, n_rows, max_levels)
      type_label <- "character"
    }

    data.frame(
      Column = col,
      Type = type_label,
      Distribution = info$html,
      Stats = info$stats,
      Complete = paste0(pct_complete, "%"),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, profiles)
}

# -- internal helpers ----------------------------------------------------------

#' @keywords internal
.profile_numeric <- function(vals, n_na, n_total) {
  clean <- vals[!is.na(vals)]
  if (length(clean) == 0) {
    return(list(
      html = "<span style='color:#999'>no data</span>",
      stats = "all NA"
    ))
  }

  # Compute stats
  q <- stats::quantile(clean, probs = c(0, 0.5, 1))
  mn <- mean(clean)

  # Build 15-bin histogram as CSS bars
  brks <- seq(q[1], q[3], length.out = 16)
  if (brks[1] == brks[16]) {
    # constant value
    html <- .css_bar_chart(rep(1, 15), "#4682B4")
  } else {
    counts <- graphics::hist(clean, breaks = brks, plot = FALSE)$counts
    html <- .css_bar_chart(counts, "#4682B4")
  }

  stats_str <- paste0(
    "min: ", round(q[1], 3),
    " | med: ", round(q[2], 3),
    " | mean: ", round(mn, 3),
    " | max: ", round(q[3], 3)
  )

  list(html = html, stats = stats_str)
}

#' @keywords internal
.profile_logical <- function(vals, n_na, n_total) {
  clean <- vals[!is.na(vals)]
  n_clean <- length(clean)
  if (n_clean == 0) {
    return(list(
      html = "<span style='color:#999'>no data</span>",
      stats = "all NA"
    ))
  }

  n_true <- sum(clean)
  n_false <- n_clean - n_true
  pct_true <- round(n_true / n_clean * 100, 1)
  pct_false <- round(100 - pct_true, 1)

  html <- paste0(
    "<div style='display:flex;height:16px;width:120px;border-radius:3px;overflow:hidden;'>",
    "<div style='width:", pct_true, "%;background:#4682B4;' title='TRUE: ",
    pct_true, "%'></div>",
    "<div style='width:", pct_false, "%;background:#D3D3D3;' title='FALSE: ",
    pct_false, "%'></div>",
    "</div>"
  )

  stats_str <- paste0("TRUE: ", n_true, " (", pct_true, "%) | FALSE: ", n_false, " (", pct_false, "%)")
  if (n_na > 0) stats_str <- paste0(stats_str, " | NA: ", n_na)

  list(html = html, stats = stats_str)
}

#' @keywords internal
.profile_categorical <- function(vals, n_total, max_levels = 5L) {
  tbl <- sort(table(vals, useNA = "no"), decreasing = TRUE)
  n_levels <- length(tbl)

  if (n_levels == 0) {
    return(list(
      html = "<span style='color:#999'>no data</span>",
      stats = "all NA"
    ))
  }

  show_levels <- utils::head(tbl, max_levels)
  palette <- c("#4682B4", "#E07941", "#50A050", "#C75D8A", "#8E6FBF")

  # Build horizontal bars for top levels
  max_count <- max(show_levels)
  bars <- vapply(seq_along(show_levels), function(i) {
    lbl <- names(show_levels)[i]
    cnt <- show_levels[i]
    pct <- round(cnt / n_total * 100, 1)
    bar_w <- round(cnt / max_count * 100)
    color <- palette[((i - 1) %% length(palette)) + 1]
    paste0(
      "<div style='display:flex;align-items:center;gap:4px;margin:1px 0;font-size:11px;line-height:1.3;'>",
      "<div style='width:", bar_w, "%;min-width:3px;height:12px;background:", color,
      ";border-radius:2px;flex-shrink:0;max-width:80px;'></div>",
      "<span style='white-space:nowrap;color:#333;'>", .html_escape(lbl),
      " (", pct, "%)</span>",
      "</div>"
    )
  }, character(1))

  html <- paste0(
    "<div style='width:200px;'>",
    paste(bars, collapse = ""),
    "</div>"
  )

  # Stats
  stats_parts <- paste0(n_levels, " unique")
  if (n_levels > max_levels) {
    stats_parts <- paste0(stats_parts, " (showing top ", max_levels, ")")
  }

  list(html = html, stats = stats_parts)
}

#' @keywords internal
.html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

#' @keywords internal
.css_bar_chart <- function(counts, color) {
  max_c <- max(counts)
  if (max_c == 0) max_c <- 1
  heights <- round(counts / max_c * 20)

  bars <- vapply(heights, function(h) {
    paste0(
      "<div style='display:inline-block;width:6px;height:", h,
      "px;background:", color, ";margin-right:1px;border-radius:1px;vertical-align:bottom;'></div>"
    )
  }, character(1))

  paste0(
    "<div style='display:inline-flex;align-items:flex-end;height:22px;'>",
    paste(bars, collapse = ""),
    "</div>"
  )
}

#' Get Seurat Object Information
#'
#' Generate a comprehensive summary of a Seurat object.
#'
#' @param seurat A Seurat object
#' @return A list with formatted information about the object
#' @export
#' @importFrom Seurat DefaultAssay VariableFeatures Idents
SeuratInfo <- function(seurat) {
  output <- list()

  # Seurat version
  output$version <- paste("Seurat version:", as.character(seurat@version))

  # Graphs
  if (length(seurat@graphs) > 0) {
    output$graphs <- paste("Graphs:", paste(names(seurat@graphs), collapse = ", "))
  } else {
    output$graphs <- "Graphs: None"
  }

  # Reductions
  if (length(seurat@reductions) > 0) {
    assayused <- sapply(names(seurat@reductions), function(name) {
      return(seurat[[name]]@assay.used)
    })
    output$reductions <- paste("Reductions:",
                               paste(names(assayused), " (", assayused, ")", sep = "", collapse = ", "))
  } else {
    output$reductions <- "Reductions: None"
  }

  # Images
  if (length(seurat@images) > 0) {
    output$images <- paste("Images:", paste(names(seurat@images), collapse = ", "))
  } else {
    output$images <- "Images: None"
  }

  # Idents
  output$ident_label <- paste("Ident label:", FindIdentLabel(seurat))

  tab <- table(Idents(seurat))
  df <- data.frame(Count = as.integer(tab))
  rownames(df) <- rownames(tab)
  output$idents_table <- t(df)

  # Assays
  assays <- names(seurat@assays)
  default_assay <- DefaultAssay(seurat)

  slotinfo <- list()
  slots <- c("counts", "data", "scale.data")

  for (assay in assays) {
    assay_obj <- seurat@assays[[assay]]

    defaultassay <- ifelse(default_assay == assay, "YES", "")

    dims <- sapply(slots, function(slot) {
      d <- get_layer_dim(assay_obj, slot)
      paste0(d[1], "x", d[2])
    })

    hvgs <- length(VariableFeatures(seurat, assay = assay))

    slotinfo[[assay]] <- c(defaultassay, dims, hvgs)
  }

  output$assays_table <- data.frame(do.call(rbind, slotinfo))
  colnames(output$assays_table) <- c("default", slots, "HVGs")

  return(output)
}
