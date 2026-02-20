#!/usr/bin/env Rscript

# Wrapper that delegates to generate_inputs.R.

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) stop("Unable to locate script path (missing --file= argument).")
  path <- sub("^--file=", "", file_arg[[1]])
  gsub("~\\+~", " ", path)
}

script_dir <- dirname(normalizePath(get_script_path()))
target <- file.path(script_dir, "generate_inputs.R")

args <- commandArgs(trailingOnly = TRUE)
status <- system2("Rscript", c(target, args))
if (!identical(status, 0L)) stop("generate_inputs.R failed with exit code: ", status)
