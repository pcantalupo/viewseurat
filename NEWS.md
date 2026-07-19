# viewseurat 0.6.2

* Add a Download TSV button to the metadata distribution plot. It exports the
  x-axis values and counts: histogram bin left edges and counts for numeric
  columns, category and count for categorical columns.

# viewseurat 0.6.1

* Add command-line scripts in `inst/scripts/`: `viewseurat.R` launches the app
  on a saved `.rds`/`.qs2` object, and `viewseurat_servers.R` lists running
  ViewSeurat servers (PID, port, object) with `--all` and `--killall` options.
* Document the command-line scripts in the README.

# viewseurat 0.6.0

* Fix viewport jump when paginating Column Profiler and Cell Metadata DT
  tables. Row-height variance across pages caused the wrapper to shrink
  dramatically on page change; a secondary browser scroll-back fired ~50 ms
  after the layout settled. Anchor the wrapper's viewport-relative bottom in
  `preDrawCallback` and restore it via double `requestAnimationFrame` plus a
  `setTimeout(100)` catch for the late scroll (#19).

# viewseurat 0.5.3

* Fix missing `metadata_utils.R` source in `app.R` that caused Column Profiler
  to fail when running via RStudio Run App.
* Add regression test guarding `app.R` source completeness.

# viewseurat 0.5.2

* Fix metadata profiler bars: scale relative to total cell count rather than
  the maximum per-category count so bar widths are proportional across
  columns (#13).

# viewseurat 0.5.1

* Fix `addResourcePath` error on startup by moving custom CSS to `inst/www`
  and loading via `system.file()`.

# viewseurat 0.5.0

* Add browser back/forward navigation support for tab and assay state.
* Add collapsible metadata column names panel to the Metadata tab.
* Reduce Docker image size by switching to `rocker/r-ver` base image.

# viewseurat 0.4.0

* Add gene search to assay tables and cell barcode search to the metadata
  table.

# viewseurat 0.3.0

* Initial tagged release.
