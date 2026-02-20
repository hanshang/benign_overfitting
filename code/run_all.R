#!/usr/bin/env Rscript

# Regenerates CSV inputs and all PNG figures referenced by main.tex.
#
# Usage:
#   Rscript code/run_all.R [--clean] [--paper] [--quick] [--force] [--reps=N] [--seed=N]

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) stop("Unable to locate script path (missing --file= argument).")
  path <- sub("^--file=", "", file_arg[[1]])
  gsub("~\\+~", " ", path)
}

script_path <- normalizePath(get_script_path())
repo_root <- normalizePath(file.path(dirname(script_path), ".."))

args <- commandArgs(trailingOnly = TRUE)

inputs_script <- file.path(repo_root, "code", "scripts", "generate_inputs.R")
figures_script <- file.path(repo_root, "code", "scripts", "generate_figures.R")

if (!file.exists(inputs_script)) stop("Missing inputs generator: ", inputs_script)
if (!file.exists(figures_script)) stop("Missing figure generator: ", figures_script)

has_flag <- function(flag) any(args == flag)

clean <- has_flag("--clean")
clean_inputs <- clean || has_flag("--clean-inputs")
clean_plots <- clean || has_flag("--clean-plots")

if (clean_inputs) {
  inputs_dir <- file.path(repo_root, "code", "inputs")
  if (dir.exists(inputs_dir)) {
    keep <- c("arcene_random.csv", "arcene_PCA.csv", "ARCENE_combined_methods.csv")
    files <- list.files(inputs_dir, full.names = TRUE, recursive = FALSE)
    to_remove <- files[!basename(files) %in% keep]
    if (length(to_remove) > 0) {
      message("Cleaning inputs (keeping ARCENE CSVs)")
      ok <- file.remove(to_remove)
      if (any(!ok)) warning("Failed to remove some input files:\n- ", paste(to_remove[!ok], collapse = "\n- "))
    }
  }
}

if (clean_plots) {
  plots_dir <- file.path(repo_root, "plots")
  figs <- c(
    "scaled_maha_gamma.png",
    "theorem1_identity.png",
    "theorem1_compoundsymmetric.png",
    "sim1.png",
    "empirical_error_identity_delta.png",
    "empirical_error_compound_symmetric.png",
    "compound_symmetric_theoretical_FINAL.png",
    "ar1_FINAL.png",
    "ar1_covariance_matrix_theoretical_FINAL.png",
    "uniform_distribution.png",
    "poisson_new_correct.png",
    "t_dist_identity_DF.png",
    "t_dist_cs_DF.png",
    "flipping_y.png",
    "signaltonoise_2.png",
    "redundant_features.png",
    "mahalanobis_distance_realdata.png",
    "arcene_random_sample.png",
    "ARCENE_Comparison.png",
    "ARCENE_PCA.png",
    "wang_maha.png",
    "wang_identity_delta.png",
    "wang_compoundsymmetric.png",
    "wang_ar1.png"
  )
  paths <- file.path(plots_dir, figs)
  existing <- paths[file.exists(paths)]
  if (length(existing) > 0) {
    ok <- file.remove(existing)
    if (any(!ok)) warning("Failed to remove some plot files:\n- ", paste(existing[!ok], collapse = "\n- "))
  }
}

arcene_script <- file.path(repo_root, "code", "scripts", "generate_arcene_inputs.R")
arcene_csvs <- file.path(
  repo_root, "code", "inputs",
  c("arcene_random.csv", "arcene_PCA.csv", "ARCENE_combined_methods.csv")
)
if (any(!file.exists(arcene_csvs))) {
  if (!file.exists(arcene_script)) stop("Missing ARCENE generator: ", arcene_script)
  message("Generating ARCENE inputs (one-time)")
  status_arcene <- system2("Rscript", arcene_script)
  if (!identical(status_arcene, 0L)) stop("ARCENE input generation failed with exit code: ", status_arcene)
}

message("Generating inputs")
status_inputs <- system2("Rscript", c(inputs_script, args))
if (!identical(status_inputs, 0L)) stop("Input generation failed with exit code: ", status_inputs)

message("\nGenerating figures")
status_figs <- system2("Rscript", c(figures_script, args))
if (!identical(status_figs, 0L)) stop("Figure generation failed with exit code: ", status_figs)
