#' Filter metadata columns by partial name match
#'
#' @param meta_data A data.frame (typically `obj@meta.data`).
#' @param search_term Character string to match against column names.
#'   Empty string or NULL returns all columns.
#' @return A data.frame with only matching columns (row names preserved).
#' @keywords internal
filter_metadata_columns <- function(meta_data, search_term) {
  search_term <- trimws(if (is.null(search_term)) "" else search_term)
  if (nchar(search_term) == 0) return(meta_data)

  matched <- grep(tolower(search_term), tolower(colnames(meta_data)), fixed = TRUE)
  if (length(matched) == 0) {
    return(meta_data[, integer(0), drop = FALSE])
  }
  meta_data[, matched, drop = FALSE]
}
