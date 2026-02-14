# CLAUDE.md

Read-only Shiny app viewer for Seurat v4 and v5 objects.

## Core Constraints
- Do NOT add analysis functions, data transformations, or object modifications.
- Use `LayerData(obj, layer = "counts")` for data access, not direct slot access

## SCTAssay Gotcha
`SCTAssay` does NOT have a `@layers` slot — direct `@layers[["..."]]` access crashes on it.
Use `get_layer_dim(assay_obj, layer)` in `seurat_utils.R` to safely get per-layer dimensions;
it checks `slotNames` first and falls back appropriately. Same pattern applies anywhere you
need to touch a layer's underlying matrix on an arbitrary assay type.

## Critical: Namespace Prefixes
ALWAYS use explicit namespace prefixes (`SeuratObject::Layers()`, `shiny::renderUI()`, `DT::datatable()`). The `app.R` launcher sources R/ files directly, bypassing NAMESPACE - code must work in that context.

## Sundry
Keep `settings.local.json` in `.claude/` directory (not project root).
Always save plan files to the project-local ./claude/plans/ directory, never to ~/.claude/plans/
Don't commit changes automatically. Wait for explicit /ship or commit request.
When implementing a plan, create a new branch.




