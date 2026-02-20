#!/usr/bin/env Rscript

# Generate all PNG figures for main.tex from CSV inputs in code/paper/inputs/.
# Output: plots/ at repo root.
# Usage: Rscript code/paper/scripts/generate_figures.R

suppressPackageStartupMessages({
  library(MASS)
})

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    stop("Unable to locate script path (missing --file= argument).")
  }
  path <- sub("^--file=", "", file_arg[[1]])
  gsub("~\\+~", " ", path)
}

script_dir <- dirname(normalizePath(get_script_path()))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
code_root <- file.path(repo_root, "code")
inputs_dir <- file.path(code_root, "inputs")
data_dir <- file.path(code_root, "data")
arcene_dir <- file.path(code_root, "ARCENE")

plots_dir <- file.path(repo_root, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  utils::read.csv(path, check.names = FALSE)
}

read_input_csv <- function(name) {
  read_csv(file.path(inputs_dir, name))
}

png_type <- if (capabilities("cairo")) "cairo" else "Xlib"

with_png <- function(path, width, height, code) {
  grDevices::png(
    filename = path,
    width = width,
    height = height,
    units = "px",
    res = 96,
    type = png_type
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  code()
  invisible(path)
}

# Wang (2018) error formula (equal class sizes), gamma = p/n in (0,1].
wang_error_equal_n <- function(delta_sq, gamma) {
  stats::pnorm(-(sqrt(1 - gamma) * (delta_sq / (2 * sqrt(delta_sq + 4 * gamma)))))
}

# Theorem 1 visualization helper. maha_scaled = Delta^2/n, gamma = p/n.
theorem1_viz <- function(maha_scaled, gamma) {
  stats::pnorm(-(2 * sqrt(maha_scaled) * sqrt(gamma * (1 - gamma))))
}

# scaled_maha_gamma.png
plot_scaled_maha_gamma <- function(out_path) {
  n <- 100
  p <- 1:100
  gamma <- p / n

  delta_id <- sqrt(p / n)

  rho_cs <- 0.1
  delta_cs <- sqrt((p / (1 - rho_cs + rho_cs * p)) / n)

  rho_ar1 <- 0.3
  one_sig_inv_one <- (p * (1 - rho_ar1) + 2 * rho_ar1) / (1 + rho_ar1)
  delta_ar1 <- sqrt(one_sig_inv_one / n)

  with_png(out_path, 800, 817, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      gamma, delta_id,
      xlab = expression(gamma),
      ylab = expression(Delta["*"]),
      main = expression(paste("Scaled Mahalanobis Distance against ", gamma)),
      pch = 19, col = "brown", type = "o", lwd = 1.2,
      xlim = c(0, 1), ylim = c(0, 1)
    )
    lines(gamma, delta_cs, col = "darkgreen", lwd = 1.2)
    points(gamma, delta_cs, pch = 19, col = "darkgreen")
    lines(gamma, delta_ar1, col = "purple", lwd = 1.2)
    points(gamma, delta_ar1, pch = 19, col = "purple")
    legend(
      x = 0.05, y = 0.95, cex = 1,
      legend = c("Identity", "Compound Symmetric", "AR(1)"),
      title = "Covariance Matrix Structure",
      col = c("brown", "darkgreen", "purple"),
      lty = 1, lwd = 1.2, bty = "n"
    )
  })
}

# theorem1_identity.png
plot_theorem1_identity <- function(out_path) {
  nk_seq <- floor(seq(20, 400, length.out = 20))
  gamma_seq <- seq(0.001, 1, length.out = 20)
  p_seq <- ceiling(nk_seq * gamma_seq)
  gamma <- p_seq / nk_seq

  deltas <- c(0.1, 0.5, 3)
  cols <- c("darkred", "black", "darkblue")
  pchs <- c(18, 16, 17)

  errs <- lapply(deltas, function(d) {
    maha_scaled <- (d^2) * (p_seq / nk_seq)
    theorem1_viz(maha_scaled, gamma)
  })

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      gamma, errs[[2]],
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Identity Covariance Matrix",
      type = "l", lwd = 1.2, col = cols[[2]],
      xlim = c(0, 1), ylim = c(0, 0.5)
    )
    for (i in seq_along(deltas)) {
      lines(gamma, errs[[i]], col = cols[[i]], lwd = 1.2)
      points(gamma, errs[[i]], col = cols[[i]], pch = pchs[[i]])
    }
    legend(
      x = 0.5, y = 0.15, cex = 1,
      legend = expression(
        paste(delta[p], " = 0.1"),
        paste(delta[p], " = 0.5"),
        paste(delta[p], " = 3")
      ),
      col = cols, lty = 1, lwd = 1.2, pch = pchs, bty = "n"
    )
  })
}

# theorem1_compoundsymmetric.png
plot_theorem1_compoundsymmetric <- function(out_path) {
  nk_seq <- floor(seq(50, 400, length.out = 20))
  gamma_seq <- seq(0.001, 1, length.out = 20)
  p_seq <- ceiling(nk_seq * gamma_seq)
  gamma <- p_seq / nk_seq

  delta <- 3
  rhos <- c(0.1, 0.3, 0.5)
  cols <- c("darkred", "black", "darkblue")
  pchs <- c(18, 16, 17)

  errs <- lapply(rhos, function(rho) {
    one_sig_inv_one <- p_seq / (1 - rho + rho * p_seq)
    maha_scaled <- (delta^2) * one_sig_inv_one / nk_seq
    theorem1_viz(maha_scaled, gamma)
  })

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      gamma, errs[[1]],
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Compound Symmetric Covariance Matrix",
      type = "l", lwd = 1.2, col = cols[[1]],
      xlim = c(0, 1), ylim = c(0, 0.5)
    )
    for (i in seq_along(rhos)) {
      lines(gamma, errs[[i]], col = cols[[i]], lwd = 1.2)
      points(gamma, errs[[i]], col = cols[[i]], pch = pchs[[i]])
    }
    legend(
      x = 0.6, y = 0.1, cex = 1,
      legend = expression(
        paste(cov(x[i], x[j]), " = 0.1"),
        paste(cov(x[i], x[j]), " = 0.3"),
        paste(cov(x[i], x[j]), " = 0.5")
      ),
      col = cols, lty = 1, lwd = 1.2, pch = pchs, bty = "n"
    )
  })
}

# sim1.png (input: sim1_n98.csv, sim1_n140.csv, sim1_n238.csv)
plot_sim1 <- function(out_path) {
  d98  <- read_input_csv("sim1_n98.csv")
  d140 <- read_input_csv("sim1_n140.csv")
  d238 <- read_input_csv("sim1_n238.csv")

  with_png(out_path, 1509, 950, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      d98$index, d98$test,
      type = "b", pch = 19, cex = 0.6,
      col = "orange", lty = 3, lwd = 1,
      xlab = "Dimension (p)", ylab = "Testing Error",
      main = "Dimension vs Testing Error (Isotropic feature)",
      xlim = c(0, 500), ylim = c(0.08, 0.45)
    )
    lines(d140$index, d140$test, type = "b", pch = 19, cex = 0.6, col = "blue", lty = 3, lwd = 1)
    lines(d238$index, d238$test, type = "b", pch = 19, cex = 0.6, col = "darkgreen", lty = 3, lwd = 1)

    abline(v = 98, col = "orange", lwd = 1)
    abline(v = 140, col = "blue", lwd = 1)
    abline(v = 238, col = "darkgreen", lwd = 1)

    legend(
      "topright", bty = "n", cex = 0.9, title = "Training Sample Size",
      legend = c("n = 98", "n = 140", "n = 238"),
      col = c("orange", "blue", "darkgreen"), lty = 3, lwd = 1
    )
  })
}

# empirical_error_identity_delta.png
plot_empirical_identity_delta <- function(out_path) {
  df <- read_input_csv("identity_empirical_gamma0to5.csv")
  del05 <- df[, 1:3]
  del1 <- df[, 4:6]
  del15 <- df[, 7:9]

  with_png(out_path, 700, 700, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    idx <- 2:nrow(del05) 
    plot(
      del05[idx, 3], del05[idx, 2],
      col = "black", type = "l",
      ylim = c(0, 0.4), xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Identity Covariance Matrix"
    )
    lines(del1[idx, 3], del1[idx, 2], col = "darkred", lwd = 1.2)
    lines(del15[idx, 3], del15[idx, 2], col = "darkblue", lwd = 1.2)

    points(del05[idx, 3], del05[idx, 2], pch = 3, col = "black")
    points(del1[idx, 3], del1[idx, 2], pch = 8, col = "darkred")
    points(del15[idx, 3], del15[idx, 2], pch = 6, col = "darkblue")

    legend(
      x = 2.5, y = 0.4, cex = 1.2,
      legend = expression(paste(delta, " = 0.5"), paste(delta, " = 1"), paste(delta, " = 1.5")),
      col = c("black", "darkred", "darkblue"),
      lty = 1, lwd = 1.2,
      pch = c(3, 8, 6),
      bty = "n"
    )
  })
}

# empirical_error_compound_symmetric.png
plot_empirical_compound_symmetric <- function(out_path) {
  df <- read_input_csv("compoundsymmetric_empirical_new.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 3], df[, 2],
      col = "#c44601", type = "l",
      ylim = c(0, 0.55), xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Compound Symmetric Covariance Matrix"
    )
    lines(df[, 3], df[, 5], col = "#5ba300", lwd = 1.2)
    lines(df[, 3], df[, 8], col = "#054fb9", lwd = 1.2)

    points(df[, 3], df[, 2], pch = 17, col = "#c44601")
    points(df[, 3], df[, 5], pch = 15, col = "#5ba300")
    points(df[, 3], df[, 8], pch = 19, col = "#054fb9")

    legend(
      x = 2, y = 0.5, cex = 1.2,
      legend = expression(
        paste(cov(x[i], x[j]), " = 0.1"),
        paste(cov(x[i], x[j]), " = 0.3"),
        paste(cov(x[i], x[j]), " = 0.5")
      ),
      col = c("#c44601", "#5ba300", "#054fb9"),
      lty = 1, lwd = 1.2, pch = c(17, 15, 19),
      bty = "n"
    )
  })
}

# compound_symmetric_theoretical_FINAL.png
plot_compound_symmetric_theoretical <- function(out_path) {
  df <- read_input_csv("theoretical_compoundsymmetric.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 1], df[, 2],
      col = "#c44601", type = "l",
      ylim = c(0, 0.55), xlim = c(0, 3), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Compound Symmetric Covariance Matrix"
    )
    lines(df[, 3], df[, 4], col = "#5ba300", lwd = 1.2)
    lines(df[, 5], df[, 6], col = "#054fb9", lwd = 1.2)

    points(df[, 1], df[, 2], pch = 17, col = "#c44601")
    points(df[, 3], df[, 4], pch = 15, col = "#5ba300")
    points(df[, 5], df[, 6], pch = 19, col = "#054fb9")

    legend(
      x = 1.8, y = 0.5, cex = 1.2,
      legend = expression(
        paste(cov(x[i], x[j]), " = 0.1"),
        paste(cov(x[i], x[j]), " = 0.3"),
        paste(cov(x[i], x[j]), " = 0.5")
      ),
      col = c("#c44601", "#5ba300", "#054fb9"),
      lty = 1, lwd = 1.2, pch = c(17, 15, 19),
      bty = "n"
    )
  })
}

# ar1_FINAL.png
plot_ar1_empirical <- function(out_path) {
  df <- read_input_csv("ar1_empirical_new.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 1], df[, 2],
      col = "#8e0bca", type = "l",
      ylim = c(0, 0.4), xlim = c(0, 3), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "AR(1) Covariance Matrix"
    )
    lines(df[, 3], df[, 4], col = "#e1111e", lwd = 1.2)
    lines(df[, 5], df[, 6], col = "#08b128", lwd = 1.2)

    points(df[, 1], df[, 2], pch = 17, col = "#8e0bca")
    points(df[, 3], df[, 4], pch = 15, col = "#e1111e")
    points(df[, 5], df[, 6], pch = 19, col = "#08b128")

    legend(
      x = 2, y = 0.4, cex = 1,
      legend = expression(paste(rho, " = 0.3"), paste(rho, " = 0.5"), paste(rho, " = 0.7")),
      col = c("#8e0bca", "#e1111e", "#08b128"),
      lty = 1, lwd = 1.2, pch = c(17, 15, 19),
      bty = "n"
    )
  })
}

# ar1_covariance_matrix_theoretical_FINAL.png
plot_ar1_theoretical <- function(out_path) {
  df <- read_input_csv("theoretical_ar1.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 1], df[, 2],
      col = "#8e0bca", type = "l",
      ylim = c(0, 0.4), xlim = c(0, 3), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "AR(1) Covariance Matrix"
    )
    lines(df[, 3], df[, 4], col = "#e1111e", lwd = 1.2)
    lines(df[, 5], df[, 6], col = "#08b128", lwd = 1.2)

    points(df[, 1], df[, 2], pch = 17, col = "#8e0bca")
    points(df[, 3], df[, 4], pch = 15, col = "#e1111e")
    points(df[, 5], df[, 6], pch = 19, col = "#08b128")

    legend(
      x = 2, y = 0.4, cex = 1,
      legend = expression(paste(rho, " = 0.3"), paste(rho, " = 0.5"), paste(rho, " = 0.7")),
      col = c("#8e0bca", "#e1111e", "#08b128"),
      lty = 1, lwd = 1.2, pch = c(17, 15, 19),
      bty = "n"
    )
  })
}

# uniform_distribution.png
plot_uniform_distribution <- function(out_path) {
  required <- file.path(inputs_dir, c("unif_05_025075.csv", "unif_007.csv", "unif_0106.csv"))
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop(
      "Missing CSV input(s) for uniform_distribution.png:\n- ",
      paste(missing, collapse = "\n- ")
    )
  }

  u1 <- read_csv(required[[1]])
  u2 <- read_csv(required[[2]])
  u3 <- read_csv(required[[3]])

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      u1[, 3], u1[, 1],
      ylim = c(0, 0.5), xlim = c(0, 4),
      main = "Uniform Distribution",
      xlab = expression(gamma), ylab = "Error Rate",
      pch = 19, col = "darkred"
    )
    lines(u1[, 3], u1[, 1], lwd = 2, col = "darkred")
    points(u2[, 3], u2[, 1], pch = 22, col = "darkgreen")
    lines(u2[, 3], u2[, 1], lwd = 2, col = "darkgreen")
    points(u3[, 3], u3[, 1], pch = 23, col = "darkblue")
    lines(u3[, 3], u3[, 1], lwd = 2, col = "darkblue")
    legend(
      x = 2, y = 0.45, bty = "n", cex = 0.8, text.col = "black",
      legend = c(
        "group 1 ~ Unif(0, 0.5), group 2 ~ Unif(0.25, 0.75)",
        "group 1 ~ Unif(0, 0.5), group 2 ~ Unif(0, 0.7)",
        "group 1 ~ Unif(0, 0.5), group 2 ~ Unif(0.1, 0.6)"
      ),
      col = c("darkred", "darkgreen", "darkblue"),
      lty = 1, pch = c(19, 22, 23)
    )
  })
}

# poisson_new_correct.png
plot_poisson_distribution <- function(out_path) {
  required <- file.path(inputs_dir, c("poisson_unif_23_34.csv", "poisson_unif_24_46.csv", "poisson_unif_67_78.csv"))
  missing <- required[!file.exists(required)]
  if (length(missing) > 0) {
    stop(
      "Missing CSV input(s) for poisson_new_correct.png:\n- ",
      paste(missing, collapse = "\n- ")
    )
  }

  p1 <- read_csv(required[[1]])
  p2 <- read_csv(required[[2]])
  p3 <- read_csv(required[[3]])

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    idx <- 1:20
    plot(
      p1[idx, 3], p1[idx, 1],
      ylim = c(0, 0.5), xlim = c(0, 4),
      main = "Poisson Distribution",
      xlab = expression(gamma), ylab = "Error Rate",
      pch = 19, col = "darkred"
    )
    lines(p1[idx, 3], p1[idx, 1], lwd = 2, col = "darkred")
    points(p2[idx, 3], p2[idx, 1], pch = 22, col = "darkgreen")
    lines(p2[idx, 3], p2[idx, 1], lwd = 2, col = "darkgreen")
    points(p3[idx, 3], p3[idx, 1], pch = 23, col = "darkblue")
    lines(p3[idx, 3], p3[idx, 1], lwd = 2, col = "darkblue")
    legend(
      x = 2, y = 0.45, bty = "n", cex = 0.8, text.col = "black",
      legend = c(
        expression(lambda[1] ~ " ~ Unif(6, 7), " ~ lambda[2] ~ " ~ Unif(7, 8)"),
        expression(lambda[1] ~ " ~ Unif(2, 4), " ~ lambda[2] ~ " ~ Unif(4, 6)"),
        expression(lambda[1] ~ " ~ Unif(2, 3), " ~ lambda[2] ~ " ~ Unif(3, 4)")
      ),
      col = c("darkblue", "darkgreen", "darkred"),
      lty = 1, pch = c(23, 22, 19)
    )
  })
}

# t_dist_identity_DF.png
plot_tdist_identity_df <- function(out_path) {
  df1 <- read_input_csv("tdist_id_df1.csv")
  df4 <- read_input_csv("tdist_id_df4.csv")
  df30 <- read_input_csv("tdist_id_df30.csv")

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df1[["gamma_result"]], df1[["result_sim1"]],
      col = "darkred", ylim = c(0, 0.5), type = "l", xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "T-Distribution Identity Covariance Matrix"
    )
    lines(df4[["gamma_result"]], df4[["result_sim1"]], col = "darkgreen", lwd = 1.2)
    lines(df30[["gamma_result"]], df30[["result_sim1"]], col = "darkblue", lwd = 1.2)
    points(df1[["gamma_result"]], df1[["result_sim1"]], pch = 17, col = "darkred")
    points(df4[["gamma_result"]], df4[["result_sim1"]], pch = 15, col = "darkgreen")
    points(df30[["gamma_result"]], df30[["result_sim1"]], pch = 19, col = "darkblue")
    legend(
      x = 2.5, y = 0.5, cex = 1.2,
      legend = c(" df = 1", " df = 4", " df = 30"),
      col = c("darkred", "darkgreen", "darkblue"),
      lty = 1, lwd = 1.2, pch = c(17, 15, 19),
      bty = "n"
    )
  })
}

# t_dist_cs_DF.png
plot_tdist_cs_df <- function(out_path) {
  df1 <- read_input_csv("cs_df1_03.csv")
  df4 <- read_input_csv("cs_df4_03.csv")
  df30 <- read_input_csv("cs_df30_03.csv")

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df1[, 3], df1[, 1],
      col = "darkred", ylim = c(0.2, 0.55), type = "l", xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "T-Distribution Compound Symmetric Covariance Matrix"
    )
    lines(df4[, 3], df4[, 1], col = "darkgreen", lwd = 1.2)
    lines(df30[, 3], df30[, 1], col = "darkblue", lwd = 1.2)
    points(df1[, 3], df1[, 1], pch = 17, col = "darkred")
    points(df4[, 3], df4[, 1], pch = 15, col = "darkgreen")
    points(df30[, 3], df30[, 1], pch = 19, col = "darkblue")
    legend(
      x = 2.5, y = 0.54, cex = 1.2,
      legend = c(" df = 1", " df = 4", " df = 30"),
      col = c("darkred", "darkgreen", "darkblue"),
      lty = 1, lwd = 1.2, pch = c(17, 15, 19),
      bty = "n"
    )
  })
}

# flipping_y.png
plot_flipping_y <- function(out_path) {
  f05 <- read_input_csv("flip_05.csv")
  f10 <- read_input_csv("flip_10.csv")
  f20 <- read_input_csv("flip_20.csv")

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      f05[, 3], f05[, 1],
      col = "darkred", ylim = c(0, 0.5), type = "p", pch = 21, xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Varying Noise"
    )
    lines(f05[, 3], f05[, 1], col = "darkred", lwd = 1.2)
    lines(f10[, 3], f10[, 1], col = "darkblue", lwd = 1.2)
    points(f10[, 3], f10[, 1], pch = 17, col = "darkblue")
    points(f20[, 3], f20[, 1], pch = 15, col = "darkgreen")
    lines(f20[, 3], f20[, 1], col = "darkgreen", lwd = 1.2)
    legend(
      x = 3, y = 0.5, title = "Noise Level",
      legend = c("5%", "10%", "20%"),
      col = c("darkred", "darkblue", "darkgreen"),
      lty = 1, lwd = 1.2, cex = 1,
      pch = c(21, 17, 15),
      bty = "n"
    )
  })
}

# signaltonoise_2.png
plot_signal_to_noise <- function(out_path) {
  s1 <- read_input_csv("signal_1_15.csv")
  s4 <- read_input_csv("signal_4_15.csv")
  s6 <- read_input_csv("signal_6_15.csv")

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      s1[, 3], s1[, 1],
      col = "darkred", ylim = c(0, 0.5), type = "p", pch = 21, xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Varying Noise"
    )
    lines(s1[, 3], s1[, 1], col = "darkred", lwd = 1.2)
    lines(s4[, 3], s4[, 1], col = "darkblue", lwd = 1.2)
    points(s4[, 3], s4[, 1], pch = 17, col = "darkblue")
    points(s6[, 3], s6[, 1], pch = 15, col = "darkgreen")
    lines(s6[, 3], s6[, 1], col = "darkgreen", lwd = 1.2)
    legend(
      x = 2.5, y = 0.5, title = "Noise Level",
      legend = expression(paste(sigma, " = 1"), paste(sigma, " = 4"), paste(sigma, " = 6")),
      col = c("darkred", "darkblue", "darkgreen"),
      lty = 1, lwd = 1.2, cex = 1,
      pch = c(21, 17, 15),
      bty = "n"
    )
  })
}

# redundant_features.png
plot_redundant_features <- function(out_path) {
  s10 <- read_input_csv("signal_4_10.csv")
  s20 <- read_input_csv("signal_4_20.csv")
  s30 <- read_input_csv("signal_4_30.csv")

  with_png(out_path, 750, 809, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      s10[, 3], s10[, 1],
      col = "darkred", ylim = c(0, 0.5), type = "p", pch = 21, xlim = c(0, 4), lwd = 1.2,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Redundant Features"
    )
    lines(s10[, 3], s10[, 1], col = "darkred", lwd = 1.2)
    lines(s20[, 3], s20[, 1], col = "darkblue", lwd = 1.2)
    points(s20[, 3], s20[, 1], pch = 17, col = "darkblue")
    points(s30[, 3], s30[, 1], pch = 15, col = "darkgreen")
    lines(s30[, 3], s30[, 1], col = "darkgreen", lwd = 1.2)
    legend(
      x = 2, y = 0.5, title = "Redundant Features Percentage",
      legend = c("10%", "20%", "30%"),
      col = c("darkred", "darkblue", "darkgreen"),
      lty = 1, lwd = 1.2, cex = 1,
      pch = c(21, 17, 15),
      bty = "n"
    )
  })
}

# mahalanobis_distance_realdata.png
plot_mahalanobis_realdata <- function(out_path, trials = 50, seed = 1) {
  train_file <- file.path(arcene_dir, "arcene_train.data")
  label_file <- file.path(arcene_dir, "arcene_train.labels")

  if (!file.exists(train_file) || !file.exists(label_file)) {
    stop(
      "ARCENE training data not found. Missing one or both files:\n- ",
      train_file,
      "\n- ",
      label_file
    )
  }

  set.seed(seed)
  y <- scan(label_file, quiet = TRUE)
  x_raw <- scan(train_file, quiet = TRUE, what = double())
  x <- matrix(x_raw, nrow = length(y), byrow = TRUE)

  g1 <- x[y == 1, , drop = FALSE]
  g2 <- x[y == -1, , drop = FALSE]

  gamma_seq <- seq(0.001, 0.9, length.out = 20)
  p_seq <- ceiling(99 * gamma_seq)
  p_seq <- p_seq[-1]
  gamma <- p_seq / 99

  maha_trials <- matrix(NA_real_, nrow = trials, ncol = length(p_seq))

  for (s in seq_len(trials)) {
    perm <- sample(ncol(x))
    for (i in seq_along(p_seq)) {
      p <- p_seq[[i]]
      idx <- perm[seq_len(p)]

      mu1 <- colMeans(g1[, idx, drop = FALSE])
      mu2 <- colMeans(g2[, idx, drop = FALSE])
      delta <- mu1 - mu2

      sigma <- (stats::cov(g1[, idx, drop = FALSE]) + stats::cov(g2[, idx, drop = FALSE])) / 2
      inv <- MASS::ginv(sigma)
      maha_trials[s, i] <- sqrt(as.numeric(t(delta) %*% inv %*% delta))
    }
  }

  maha_mean <- colMeans(maha_trials, na.rm = TRUE)

  with_png(out_path, 772, 612, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      gamma, maha_mean,
      type = "b", pch = 1, lwd = 1.2,
      xlab = expression(gamma), ylab = "Mahalanobis Distance",
      main = "Mahalanobis Distance"
    )
  })
}

# arcene_random_sample.png
plot_arcene_random_sample <- function(out_path) {
  df <- read_input_csv("arcene_random.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[[1]], df[[2]],
      type = "b", pch = 19, cex = 0.7,
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Arcene Data (Randomly Sample Features)",
      xlim = c(0, 10)
    )
    abline(v = 1, col = "red")
  })
}

# ARCENE_Comparison.png
plot_arcene_comparison <- function(out_path) {
  df <- read_input_csv("ARCENE_combined_methods.csv")
  gamma <- df[["V1"]]

  with_png(out_path, 771, 759, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      gamma, df[["error"]],
      type = "o", pch = 19, cex = 0.7,
      col = "darkred", lwd = 1.2,
      ylim = c(0.1, 0.6),
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Comparison of Error Rate"
    )
    lines(gamma, df[["dlda"]], type = "o", pch = 19, cex = 0.7, col = "blue", lwd = 1.2)
    lines(gamma, df[["slda"]], type = "o", pch = 19, cex = 0.7, col = "orange", lwd = 1.2)
    lines(gamma, df[["mlda"]], type = "o", pch = 19, cex = 0.7, col = "darkgreen", lwd = 1.2)
    legend(
      "topright", bty = "n",
      legend = c("Pseudo-LDA", "DLDA", "SLDA", "MLDA"),
      col = c("darkred", "blue", "orange", "darkgreen"),
      lty = 1, lwd = 1.2, pch = 19
    )
  })
}

# ARCENE_PCA.png
plot_arcene_pca <- function(out_path) {
  df <- read_input_csv("arcene_PCA.csv")
  gamma <- df[[1]]
  pcs <- ceiling(gamma * 99)
  err <- df[["test_error"]]

  # Subsample every 5th point.
  idx <- seq(1, length(pcs), by = 5)
  pcs_sub <- pcs[idx]
  err_sub <- err[idx]

  with_png(out_path, 800, 800, function() {
    par(mar = c(5, 5, 2, 2) + 0.1)
    plot(
      pcs_sub, err_sub,
      pch = 16, cex = 0.7,
      ylim = c(0.1, 0.3),
      xlab = "Principal Components", ylab = "Error Rate"
    )
    lines(lowess(pcs_sub, err_sub), col = "darkblue", lwd = 1.8)
  })
}

# wang_maha.png
plot_wang_maha <- function(out_path) {
  gamma <- seq(0, 1, length.out = 100)
  deltas <- c(0.5, 2, 5)
  cols <- c("black", "blue", "darkred")
  pchs <- c(16, 17, 18)

  errs <- lapply(deltas, function(d2) wang_error_equal_n(d2, gamma))

  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      gamma, errs[[1]],
      xlab = expression(gamma), ylab = "Error Rate",
      main = "Mahalanobis Distance",
      type = "l", lwd = 1.2, col = cols[[1]],
      ylim = c(0, 0.5), xlim = c(0, 1)
    )
    for (i in seq_along(deltas)) {
      lines(gamma, errs[[i]], col = cols[[i]], lwd = 1.2)
      points(gamma[seq(1, length(gamma), by = 10)], errs[[i]][seq(1, length(gamma), by = 10)], pch = pchs[[i]], col = cols[[i]])
    }
    legend(
      x = 0.55, y = 0.2, cex = 1,
      legend = expression(paste(Delta^2, " = 0.5"), paste(Delta^2, " = 2"), paste(Delta^2, " = 5")),
      col = cols, lty = 1, lwd = 1.2, pch = pchs,
      bty = "n"
    )
  })
}

# wang_identity_delta.png
plot_wang_identity_delta <- function(out_path) {
  df <- read_input_csv("wang_identity_delta.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 1], df[, 2],
      xlab = expression(gamma), ylab = "Error Rate",
      type = "l", ylim = c(0, 0.5), lwd = 1.2,
      main = "Identity Covariance Matrix"
    )
    lines(df[, 3], df[, 4], col = "darkblue", lwd = 1.2)
    lines(df[, 5], df[, 6], col = "darkred", lwd = 1.2)
    points(df[, 1], df[, 2], pch = 16)
    points(df[, 3], df[, 4], col = "darkblue", pch = 17)
    points(df[, 5], df[, 6], col = "darkred", pch = 18)
    legend(
      x = 0.4, y = 0.5, cex = 1,
      legend = expression(paste(delta, " = 0.5"), paste(delta, " = 1"), paste(delta, " = 2")),
      col = c("black", "darkblue", "darkred"),
      lty = 1, lwd = 1.2, pch = c(16, 17, 18),
      bty = "n"
    )
  })
}

# wang_compoundsymmetric.png
plot_wang_compoundsymmetric <- function(out_path) {
  df <- read_input_csv("wang_compoundsymmetric.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 1], df[, 2],
      xlab = expression(gamma), ylab = "Error Rate",
      type = "l", ylim = c(0, 0.5), lwd = 1.2,
      main = "Compound Symmetric Covariance Matrix"
    )
    lines(df[, 3], df[, 4], col = "darkblue", lwd = 1.2)
    lines(df[, 5], df[, 6], col = "darkred", lwd = 1.2)
    points(df[, 1], df[, 2], pch = 16)
    points(df[, 3], df[, 4], col = "darkblue", pch = 17)
    points(df[, 5], df[, 6], col = "darkred", pch = 18)
    legend(
      x = 0.6, y = 0.1, cex = 1,
      legend = expression(
        paste(cov(X[i], X[j]), " = 0.1"),
        paste(cov(X[i], X[j]), " = 0.3"),
        paste(cov(X[i], X[j]), " = 0.5")
      ),
      col = c("black", "darkblue", "darkred"),
      lty = 1, lwd = 1.2, pch = c(16, 17, 18),
      bty = "n"
    )
  })
}

# wang_ar1.png
plot_wang_ar1 <- function(out_path) {
  df <- read_input_csv("wang_ar1.csv")
  with_png(out_path, 750, 750, function() {
    par(mar = c(5, 5, 4, 2) + 0.1)
    plot(
      df[, 1], df[, 2],
      xlab = expression(gamma), ylab = "Error Rate",
      type = "l", ylim = c(0, 0.5), lwd = 1.2,
      main = "AR(1) Covariance Matrix"
    )
    lines(df[, 3], df[, 4], col = "darkblue", lwd = 1.2)
    lines(df[, 5], df[, 6], col = "darkred", lwd = 1.2)
    points(df[, 1], df[, 2], pch = 16)
    points(df[, 3], df[, 4], col = "darkblue", pch = 17)
    points(df[, 5], df[, 6], col = "darkred", pch = 18)
    legend(
      x = 0.4, y = 0.5, cex = 1,
      legend = expression(
        paste(cov(x[i], x[j]), " = 0.1"),
        paste(cov(x[i], x[j]), " = 0.3"),
        paste(cov(x[i], x[j]), " = 0.5")
      ),
      col = c("black", "darkblue", "darkred"),
      lty = 1, lwd = 1.2, pch = c(16, 17, 18),
      bty = "n"
    )
  })
}

out_files <- c(
  plot_scaled_maha_gamma(file.path(plots_dir, "scaled_maha_gamma.png")),
  plot_theorem1_identity(file.path(plots_dir, "theorem1_identity.png")),
  plot_theorem1_compoundsymmetric(file.path(plots_dir, "theorem1_compoundsymmetric.png")),
  plot_sim1(file.path(plots_dir, "sim1.png")),
  plot_empirical_identity_delta(file.path(plots_dir, "empirical_error_identity_delta.png")),
  plot_empirical_compound_symmetric(file.path(plots_dir, "empirical_error_compound_symmetric.png")),
  plot_compound_symmetric_theoretical(file.path(plots_dir, "compound_symmetric_theoretical_FINAL.png")),
  plot_ar1_empirical(file.path(plots_dir, "ar1_FINAL.png")),
  plot_ar1_theoretical(file.path(plots_dir, "ar1_covariance_matrix_theoretical_FINAL.png")),
  plot_uniform_distribution(file.path(plots_dir, "uniform_distribution.png")),
  plot_poisson_distribution(file.path(plots_dir, "poisson_new_correct.png")),
  plot_tdist_identity_df(file.path(plots_dir, "t_dist_identity_DF.png")),
  plot_tdist_cs_df(file.path(plots_dir, "t_dist_cs_DF.png")),
  plot_flipping_y(file.path(plots_dir, "flipping_y.png")),
  plot_signal_to_noise(file.path(plots_dir, "signaltonoise_2.png")),
  plot_redundant_features(file.path(plots_dir, "redundant_features.png")),
  plot_mahalanobis_realdata(file.path(plots_dir, "mahalanobis_distance_realdata.png")),
  plot_arcene_random_sample(file.path(plots_dir, "arcene_random_sample.png")),
  plot_arcene_comparison(file.path(plots_dir, "ARCENE_Comparison.png")),
  plot_arcene_pca(file.path(plots_dir, "ARCENE_PCA.png")),
  plot_wang_maha(file.path(plots_dir, "wang_maha.png")),
  plot_wang_identity_delta(file.path(plots_dir, "wang_identity_delta.png")),
  plot_wang_compoundsymmetric(file.path(plots_dir, "wang_compoundsymmetric.png")),
  plot_wang_ar1(file.path(plots_dir, "wang_ar1.png"))
)

message("Wrote figures to: ", plots_dir)
invisible(unlist(out_files))
