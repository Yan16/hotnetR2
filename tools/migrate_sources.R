# Mechanical migration of the frozen source copies; run once during development.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(args[[1]], mustWork = TRUE)
files <- c(list.files(file.path(root, "R"), "[.]R$", full.names = TRUE),
           list.files(file.path(root, "tests"), "[.]R$", full.names = TRUE, recursive = TRUE),
           list.files(file.path(root, "inst", "scripts"), "[.]R$", full.names = TRUE))
for (file in files) {
  text <- paste(readLines(file, warn = FALSE), collapse = "\n")
  if (grepl("workflow_", basename(file)) || grepl("/tests/", file)) {
    text <- gsub("summarize_ldak_results", "summarize_analysis_ldak_results", text, fixed = TRUE)
    text <- gsub("filter_network_ldak(", "filter_analysis_network_ldak(", text, fixed = TRUE)
    text <- gsub("filter_network_ldak <-", "filter_analysis_network_ldak <-", text, fixed = TRUE)
  }
  text <- gsub("hotnetR:::", "", text, fixed = TRUE)
  text <- gsub("hotnetR::", "", text, fixed = TRUE)
  text <- gsub("hotnethelper", "hotnetR2", text, fixed = TRUE)
  text <- gsub('package = "hotnetR"', 'package = "hotnetR2"', text, fixed = TRUE)
  text <- gsub('library(hotnetR)', 'library(hotnetR2)', text, fixed = TRUE)
  text <- gsub("#' @import tidyverse\n", "", text, fixed = TRUE)
  writeLines(text, file)
}
# Remove the earlier duplicate; r_utils.R holds the last-loaded legacy definition.
file <- file.path(root, "R", "data.R")
lines <- readLines(file)
expr <- parse(file, keep.source = TRUE)
first_end <- attr(expr, "srcref")[[1]][[3]]
writeLines(lines[-seq_len(first_end)], file)

# Qualify tidyverse calls without changing local functions or base/stats calls.
packages <- c("dplyr", "tidyr", "readr", "purrr", "tibble", "stringr", "rlang", "tidyselect",
              "igraph", "yaml", "utils", "tools", "grDevices", "gtexr", "readxl")
function_map <- list()
for (package in rev(packages)) {
  for (name in getNamespaceExports(package)) function_map[[name]] <- package
}
base_names <- c(ls(baseenv(), all.names = TRUE), getNamespaceExports("stats"))
local_names <- character()
for (file in list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)) {
  expr <- parse(file)
  for (x in as.list(expr)) {
    if (is.call(x) && identical(x[[1]], as.name("<-")) && is.symbol(x[[2]])) {
      local_names <- c(local_names, as.character(x[[2]]))
    }
  }
}
for (file in list.files(file.path(root, "R"), "[.]R$", full.names = TRUE)) {
  text <- paste(readLines(file), collapse = "\n")
  tokens <- utils::getParseData(parse(text = text, keep.source = TRUE))
  calls <- tokens[tokens$token == "SYMBOL_FUNCTION_CALL" &
                    !tokens$text %in% c(base_names, local_names), ]
  # Source-token edits retain comments, Roxygen, argument names and formatting.
  lines <- strsplit(text, "\n", fixed = TRUE)[[1]]
  calls <- calls[order(calls$line1, calls$col1, decreasing = TRUE), ]
  for (i in seq_len(nrow(calls))) {
    token <- calls[i, ]
    package <- function_map[[token$text]]
    if (is.null(package)) next
    line <- lines[[token$line1]]
    prefix <- substr(line, 1L, token$col1 - 1L)
    if (grepl("[:$@]\\s*$", prefix)) next
    lines[[token$line1]] <- paste0(prefix, package, "::", substr(line, token$col1, nchar(line)))
  }
  writeLines(lines, file)
}
