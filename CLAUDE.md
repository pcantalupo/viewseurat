# CLAUDE.md

Read-only Shiny app viewer for Seurat v4 and v5 objects.

## Core Constraints
- Do NOT add analysis functions, data transformations, or object modifications.
- Use `LayerData(obj, layer = "counts")` for data access, not direct slot access

## Seurat/Bioinformatics
- `SCTAssay` does NOT have a `@layers` slot — direct `@layers[["..."]]` access crashes on it.
- Use `get_layer_dim(assay_obj, layer)` in `seurat_utils.R` to safely get per-layer dimensions; it checks `slotNames` first and falls back appropriately. Same pattern applies anywhere you need to touch a layer's underlying matrix on an arbitrary assay type.
- Seurat uses specialized assay classes: ChromatinAssay (Signac) for ATAC, SCTAssay for SCT
- VisiumV2 is a subclass of FOV

## Shiny/CSS Gotcha
AdminLTE has deep CSS specificity — CSS overrides must target actual shinydashboard DOM structure with browser-inspector-level selectors.

## Critical: Namespace Prefixes
ALWAYS use explicit namespace prefixes (`SeuratObject::Layers()`, `shiny::renderUI()`, `DT::datatable()`). The `app.R` launcher sources R/ files directly, bypassing NAMESPACE.

## Sundry
- Keep `settings.local.json` in `.claude/` directory (not project root).





