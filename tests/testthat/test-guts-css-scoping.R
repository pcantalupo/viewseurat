# Tests for Guts tab CSS scoping -----------------------------------------------
#
# Regression: tabsets.min.css was loaded globally in <head>, so its unscoped
# `.tab-pane` rule added a faint border to Bootstrap/AdminLTE's own .tab-pane
# elements on every page of the app.
#
# Fix: CSS lives in www/guts-scoped.css with every selector prefixed by
# .guts-tabset, and is loaded only when the Guts tab renders.

css_path <- function() test_path("../../www/guts-scoped.css")

# -- www/guts-scoped.css exists and is scoped ----------------------------------

test_that("www/guts-scoped.css exists", {
  expect_true(file.exists(css_path()))
})

test_that("guts-scoped.css: every .tab-pane rule is under .guts-tabset (regression guard)", {
  skip_if_not(file.exists(css_path()))

  lines <- readLines(css_path())
  rule_lines <- lines[!grepl("^\\s*(/\\*|\\*)", lines)]  # skip comment lines
  lines_with_tab_pane <- rule_lines[grepl("\\.tab-pane", rule_lines)]

  expect_gt(length(lines_with_tab_pane), 0L)
  expect_true(
    all(grepl("\\.guts-tabset", lines_with_tab_pane)),
    label = "bare .tab-pane found — would pollute Bootstrap .tab-pane on every app page"
  )
})

test_that("guts-scoped.css: every .tab-link rule is under .guts-tabset", {
  skip_if_not(file.exists(css_path()))

  lines <- readLines(css_path())
  rule_lines <- lines[!grepl("^\\s*(/\\*|\\*)", lines)]  # skip comment lines
  lines_with_tab_link <- rule_lines[grepl("\\.tab-link", rule_lines)]

  expect_gt(length(lines_with_tab_link), 0L)
  expect_true(
    all(grepl("\\.guts-tabset", lines_with_tab_link)),
    label = "bare .tab-link found — would interfere with Bootstrap nav tabs globally"
  )
})

test_that("guts-scoped.css: every .hljs rule is under .guts-tabset", {
  skip_if_not(file.exists(css_path()))

  lines <- readLines(css_path())
  rule_lines <- lines[!grepl("^\\s*(/\\*|\\*)", lines)]  # skip comment lines
  lines_with_hljs <- rule_lines[grepl("\\.hljs", rule_lines)]

  expect_gt(length(lines_with_hljs), 0L)
  expect_true(
    all(grepl("\\.guts-tabset", lines_with_hljs)),
    label = "bare .hljs rule found — highlight.js styles should be scoped to Guts tab"
  )
})

# -- viewseurat_ui() global head is clean --------------------------------------

test_that("viewseurat_ui() does not inject tabsets CDN stylesheet globally", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinydashboard")

  html <- as.character(viewseurat:::viewseurat_ui())

  expect_false(
    grepl("tabsets\\.min\\.css", html),
    label = "tabsets.min.css in global <head> leaks .tab-pane border onto every page"
  )
})

test_that("viewseurat_ui() does not inject highlight.js CDN stylesheet globally", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinydashboard")

  html <- as.character(viewseurat:::viewseurat_ui())

  expect_false(
    grepl("highlight\\.min\\.css", html),
    label = "highlight.min.css in global <head> leaks syntax-highlight styles globally"
  )
})

test_that("viewseurat_ui() does not inject highlight.js CDN script globally", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("shinydashboard")

  html <- as.character(viewseurat:::viewseurat_ui())

  expect_false(
    grepl("highlight\\.min\\.js", html),
    label = "highlight.min.js loaded in global <head> — should load on-demand in Guts tab only"
  )
})
