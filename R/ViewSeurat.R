#' View a Seurat Object Interactively
#'
#' Launch an interactive Shiny dashboard to explore a Seurat object.
#' When called with an object, the Upload tab is hidden and the viewer
#' navigates directly to the Overview tab. When called without an object,
#' the full upload interface is shown.
#'
#' @param obj A Seurat object to view. If NULL, launches the app with
#'   the upload interface visible.
#' @param title Character string for the application title displayed in the
#'   dashboard header. If NULL (default), no custom title is shown. When a
#'   file is uploaded via the interface, the filename is displayed automatically.
#' @param launch.browser Logical indicating whether to open the app in
#'   a browser window. Default is TRUE.
#' @param port The TCP port that the application should listen on. If NULL
#'   (the default), a random port will be chosen.
#' @param host The IPv4 address that the application should listen on.
#'   Default is "127.0.0.1" (localhost only).
#'
#' @return This function does not return a value. It launches a Shiny
#'   application in the viewer or browser.
#'
#' @examples
#' \dontrun{
#' # View a Seurat object directly
#' ViewSeurat(seurat_obj)
#'
#' # View with custom title
#' ViewSeurat(seurat_obj, title = "My Analysis")
#'
#' # Launch without an object (shows upload interface)
#' ViewSeurat()
#'
#' # Specify port
#' ViewSeurat(seurat_obj, port = 3838)
#' }
#'
#' @export
ViewSeurat <- function(obj = NULL,
                       title = NULL,
                       launch.browser = TRUE,
                       port = NULL,
                       host = "127.0.0.1") {

  # Validate Seurat object if provided
  if (!is.null(obj)) {
    if (!inherits(obj, "Seurat")) {
      stop("Object must be a Seurat object", call. = FALSE)
    }
    validate_seurat_object(obj)
  }

  # Pass object and state to the Shiny app via shinyOptions
  shiny::shinyOptions(viewseurat.obj = obj)
  shiny::shinyOptions(viewseurat.preloaded = !is.null(obj))
  shiny::shinyOptions(viewseurat.title = title)

  # Create the app object
  app <- viewseurat_app()

  # Build run arguments
  run_args <- list(
    appDir = app,
    launch.browser = launch.browser,
    host = host
  )

  if (!is.null(port)) {
    run_args$port <- port
  }

  # Launch the app
  do.call(shiny::runApp, run_args)
}

