#!/usr/bin/env Rscript

pacman::p_load('optparse')

##################### Options ########################
opt_parser <- OptionParser(
  option_list=list(
    make_option("--seurat", default = "seurat.qs2", type = "character",
                help = "Seurat object RDS or QS2 [default: %default]"),
    make_option("--title", default = "Seurat object", type = "character",
                help = "Title to display in the ViewSeurat app [default: %default]"),
    make_option("--port", default = 5243, type = "integer",
                help="Port number [default: %default]"),
    make_option("--launch", default = TRUE, type = "logical",
                help="Launch browser? [default: %default]")
  )
)
opts <- parse_args(opt_parser)

if (!file.exists(opts$seurat)) {
  print_help(opt_parser)
  stop("Seurat object ", opts$seurat, " was not found")
}
####################################################

pacman::p_load(qs2, viewseurat, Seurat)

# Detect file format and read appropriately
if (grepl("\\.qs2?$", opts$seurat, ignore.case = TRUE)) {
  seurat <- qs_read(opts$seurat)
} else {
  seurat <- readRDS(opts$seurat)
}

ViewSeurat(seurat, title = opts$title, port = opts$port, launch.browser = opts$launch)





