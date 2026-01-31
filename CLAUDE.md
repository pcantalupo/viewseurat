# CLAUDE.md

Read-only Shiny app viewer for Seurat v5 objects.

## Core Constraints
- Do NOT add analysis functions, data transformations, or object modifications.
- Seurat v5 only: Use `GetAssayData(obj, layer = "counts")` for data access, not direct slot access

## Critical: Namespace Prefixes
ALWAYS use explicit namespace prefixes (`SeuratObject::Layers()`, `shiny::renderUI()`, `DT::datatable()`). The `app.R` launcher sources R/ files directly, bypassing NAMESPACE - code must work in that context.

