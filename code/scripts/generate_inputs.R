#!/usr/bin/env Rscript

# Generate CSV inputs for generate_figures.R via Monte Carlo simulation.

options(stringsAsFactors = FALSE)

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) stop("Unable to locate script path (missing --file= argument).")
  path <- sub("^--file=", "", file_arg[[1]])
  gsub("~\\+~", " ", path)
}

script_dir <- dirname(normalizePath(get_script_path()))
repo_root <- normalizePath(file.path(script_dir, "..", ".."))
code_root <- file.path(repo_root, "code")
inputs_dir <- file.path(code_root, "inputs")
data_dir <- file.path(code_root, "data")

dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "experiments_lib.R"))

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(flag) any(args == flag)
get_arg_value <- function(prefix) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(NULL)
  sub(paste0("^", prefix), "", hit[[1]])
}

quick <- has_flag("--quick")
force <- has_flag("--force")
paper <- has_flag("--paper")

only_arg <- get_arg_value("--only=")
only_sections <- NULL
available_sections <- c(
  "sim1",
  "theory",
  "identity",
  "cs",
  "ar1",
  "uniform",
  "signal",
  "tdist_id",
  "poisson",
  "tdist_cs",
  "flip"
)
if (!is.null(only_arg)) {
  only_sections <- strsplit(only_arg, ",", fixed = TRUE)[[1]]
  only_sections <- trimws(only_sections)
  only_sections <- only_sections[nzchar(only_sections)]
  if (length(only_sections) == 0) stop("Invalid --only= value (empty).")
  unknown <- setdiff(tolower(only_sections), available_sections)
  if (length(unknown) > 0) {
    stop(
      "Unknown --only section(s): ",
      paste(unknown, collapse = ", "),
      "\nValid: ",
      paste(available_sections, collapse = ", ")
    )
  }
}

should_run <- function(section) {
  if (is.null(only_sections)) return(TRUE)
  any(tolower(trimws(section)) == tolower(only_sections))
}

seed <- 123L
seed_arg <- get_arg_value("--seed=")
if (!is.null(seed_arg)) {
  seed <- suppressWarnings(as.integer(seed_arg))
  if (is.na(seed)) stop("Invalid --seed value: ", seed_arg)
}

reps_arg <- get_arg_value("--reps=")
reps_override <- NULL
if (!is.null(reps_arg)) {
  reps_override <- suppressWarnings(as.integer(reps_arg))
  if (is.na(reps_override) || reps_override < 1) stop("Invalid --reps value: ", reps_arg)
}

# Choose MC repetition count:
resolve_reps <- function(paper_default) {
  if (!is.null(reps_override)) return(as.integer(reps_override))
  if (quick) return(8L)
  if (paper) return(as.integer(paper_default))
  40L
}

# Per-section repetition counts.
reps_identity <- resolve_reps(100L)
reps_cs <- resolve_reps(100L)
reps_ar1 <- resolve_reps(100L)
reps_unif <- resolve_reps(100L)
reps_poisson <- resolve_reps(100L)
reps_signal <- resolve_reps(80L)
reps_tdist_id <- resolve_reps(100L)
reps_tdist_cs <- resolve_reps(100L)
reps_flip <- resolve_reps(100L)
reps_sim1 <- resolve_reps(50L)

pipeline_version <- 3L

write_csv <- function(df, path) utils::write.csv(df, path, row.names = FALSE)

meta_path_for <- function(path) paste0(path, ".meta.rds")

# Check if output needs regeneration.
need_output <- function(path, meta) {
  if (force || !file.exists(path)) return(TRUE)
  mp <- meta_path_for(path)
  if (!file.exists(mp)) return(FALSE)
  old <- tryCatch(readRDS(mp), error = function(e) NULL)
  if (is.null(old)) return(TRUE)
  !isTRUE(all.equal(old, meta, check.attributes = FALSE))
}

# Total training samples (both classes) after stratified split.
train_total_for <- function(n_per_class, train_frac) {
  per_class <- ceiling(train_frac * n_per_class)
  per_class <- max(1L, min(n_per_class - 1L, per_class))
  2L * as.integer(per_class)
}

# MC simulation of LDA test error over a grid of (n, p) pairs.
simulate_lda_curve <- function(
  n_per_class,
  p_dim,
  sampler_pair,
  train_frac = 0.7,
  reps = 40,
  seed = 1L,
  fit_fn = lda_ginv_svd
) {
  stopifnot(length(n_per_class) == length(p_dim))
  if (!is.numeric(train_frac) || length(train_frac) != 1 || train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a scalar in (0, 1).")
  }
  if (!is.numeric(reps) || length(reps) != 1 || reps < 1) stop("reps must be >= 1.")

  set.seed(seed)
  out_mean <- numeric(length(n_per_class))
  out_sd <- numeric(length(n_per_class))
  out_gamma <- numeric(length(n_per_class))

  for (i in seq_along(n_per_class)) {
    n <- as.integer(n_per_class[[i]])
    p <- max(1L, as.integer(p_dim[[i]]))
    y <- c(rep(1, n), rep(-1, n))

    errs <- numeric(reps)
    train_total <- train_total_for(n, train_frac)

    for (r in seq_len(reps)) {
      samp <- sampler_pair(n, p)
      x <- rbind(samp$g1, samp$g2)

      split <- split_train_test(y, train_frac = train_frac)
      train_idx <- split$train
      test_idx <- split$test

      mod <- fit_fn(x[train_idx, , drop = FALSE], y[train_idx])
      pred <- predict_lda(x[test_idx, , drop = FALSE], mod)
      errs[[r]] <- mean(pred != y[test_idx])
    }

    out_mean[[i]] <- mean(errs)
    out_sd[[i]] <- if (length(errs) > 1) stats::sd(errs) else NA_real_
    out_gamma[[i]] <- p / train_total

    message(sprintf(
      "point %02d/%02d: n=%d p=%d gamma=%.3f mean=%.4f sd=%.4f",
      i, length(n_per_class), n, p, out_gamma[[i]], out_mean[[i]], out_sd[[i]]
    ))
  }

  data.frame(result_sim1 = out_mean, result_se = out_sd, gamma_result = out_gamma)
}

# MC simulation returning both train and test error.
simulate_lda_curve_train_test <- function(n_per_class, p_dim, sampler_pair, train_frac = 0.7, reps = 40, seed = 1L) {
  stopifnot(length(n_per_class) == length(p_dim))
  if (!is.numeric(train_frac) || length(train_frac) != 1 || train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a scalar in (0, 1).")
  }
  if (!is.numeric(reps) || length(reps) != 1 || reps < 1) stop("reps must be >= 1.")

  set.seed(seed)
  out_train <- numeric(length(n_per_class))
  out_test <- numeric(length(n_per_class))
  out_gamma <- numeric(length(n_per_class))

  for (i in seq_along(n_per_class)) {
    n <- as.integer(n_per_class[[i]])
    p <- max(1L, as.integer(p_dim[[i]]))
    y <- c(rep(1, n), rep(-1, n))
    train_total <- train_total_for(n, train_frac)

    train_errs <- numeric(reps)
    test_errs <- numeric(reps)

    for (r in seq_len(reps)) {
      samp <- sampler_pair(n, p)
      x <- rbind(samp$g1, samp$g2)

      split <- split_train_test(y, train_frac = train_frac)
      train_idx <- split$train
      test_idx <- split$test

      mod <- lda_ginv_svd(x[train_idx, , drop = FALSE], y[train_idx])
      pred_train <- predict_lda(x[train_idx, , drop = FALSE], mod)
      pred_test <- predict_lda(x[test_idx, , drop = FALSE], mod)
      train_errs[[r]] <- mean(pred_train != y[train_idx])
      test_errs[[r]] <- mean(pred_test != y[test_idx])
    }

    out_train[[i]] <- mean(train_errs)
    out_test[[i]] <- mean(test_errs)
    out_gamma[[i]] <- p / train_total
    message(sprintf(
      "point %02d/%02d: n=%d p=%d gamma=%.3f train=%.4f test=%.4f",
      i, length(n_per_class), n, p, out_gamma[[i]], out_train[[i]], out_test[[i]]
    ))
  }

  data.frame(train_error = out_train, test_error = out_test, gamma_result = out_gamma)
}

# MC simulation with random label flipping (identity cov).
simulate_lda_curve_label_flip <- function(n_per_class, p_dim, delta, flip_frac, train_frac = 0.7, reps = 40, seed = 1L) {
  stopifnot(length(n_per_class) == length(p_dim))
  if (!is.numeric(delta) || length(delta) != 1) stop("delta must be a scalar.")
  if (!is.numeric(flip_frac) || length(flip_frac) != 1 || flip_frac < 0 || flip_frac > 0.5) {
    stop("flip_frac must be a scalar in [0, 0.5].")
  }

  set.seed(seed)
  out_mean <- numeric(length(n_per_class))
  out_sd <- numeric(length(n_per_class))
  out_gamma <- numeric(length(n_per_class))

  for (i in seq_along(n_per_class)) {
    n <- as.integer(n_per_class[[i]])
    p <- max(1L, as.integer(p_dim[[i]]))
    train_total <- train_total_for(n, train_frac)

    errs <- numeric(reps)
    for (r in seq_len(reps)) {
      mu1 <- rep(0, p)
      mu2 <- rep(delta, p)
      g1 <- sample_normal_identity(n, mu1)
      g2 <- sample_normal_identity(n, mu2)
      x <- rbind(g1, g2)
      y <- c(rep(1, n), rep(-1, n))

      n_flip <- max(0L, as.integer(ceiling(flip_frac * length(y))))
      if (n_flip > 0) {
        idx <- sample.int(length(y), size = n_flip, replace = FALSE)
        y[idx] <- -y[idx]
      }

      split <- split_train_test(y, train_frac = train_frac)
      train_idx <- split$train
      test_idx <- split$test

      mod <- lda_ginv_svd(x[train_idx, , drop = FALSE], y[train_idx])
      pred <- predict_lda(x[test_idx, , drop = FALSE], mod)
      errs[[r]] <- mean(pred != y[test_idx])
    }

    out_mean[[i]] <- mean(errs)
    out_sd[[i]] <- if (length(errs) > 1) stats::sd(errs) else NA_real_
    out_gamma[[i]] <- p / train_total
    message(sprintf(
      "point %02d/%02d: n=%d p=%d gamma=%.3f flip=%.0f%% mean=%.4f sd=%.4f",
      i, length(n_per_class), n, p, out_gamma[[i]], 100 * flip_frac, out_mean[[i]], out_sd[[i]]
    ))
  }

  data.frame(result_sim1 = out_mean, result_se = out_sd, gamma_result = out_gamma)
}

# sim1: dimension vs testing error

if (should_run("sim1")) {
  sim1_p_seq <- seq(1L, 481L, by = 5L)
  sim1_configs <- list(
    n98  = list(n_per_class = 70L,  label = "n98"),
    n140 = list(n_per_class = 100L, label = "n140"),
    n238 = list(n_per_class = 170L, label = "n238")
  )
  sim1_delta <- 0.3

  for (cfg_name in names(sim1_configs)) {
    cfg <- sim1_configs[[cfg_name]]
    out_path <- file.path(inputs_dir, paste0("sim1_", cfg_name, ".csv"))

    meta <- list(
      version = pipeline_version,
      kind = "sim1_dimension_vs_error",
      delta = sim1_delta,
      n_per_class = cfg$n_per_class,
      p_seq = sim1_p_seq,
      train_frac = 0.7,
      reps = reps_sim1,
      seed = as.integer(seed + 500L)
    )

    if (!need_output(out_path, meta)) next
    message("\nGenerating sim1 ", cfg_name, " (Monte Carlo, n_per_class=", cfg$n_per_class, ")")

    n_pc <- cfg$n_per_class
    res <- simulate_lda_curve_train_test(
      n_per_class = rep(n_pc, length(sim1_p_seq)),
      p_dim = sim1_p_seq,
      sampler_pair = function(n, p) {
        mu1 <- rep(0, p)
        mu2 <- rep(sim1_delta, p)
        list(
          g1 = sample_normal_identity(n, mu1),
          g2 = sample_normal_identity(n, mu2)
        )
      },
      train_frac = 0.7,
      reps = reps_sim1,
      seed = seed + 500L + n_pc
    )

    out_df <- data.frame(
      index = sim1_p_seq,
      train = res$train_error,
      test  = res$test_error
    )
    write_csv(out_df, out_path)
    saveRDS(meta, meta_path_for(out_path))
    message("Wrote: ", out_path)
  }
}

reps_theory <- resolve_reps(50L)

# Mahalanobis distance 1'Sigma^{-1}1 for AR(1) with Sigma_ij = rho^|i-j|.
ar1_ones_maha <- function(p, rho) {
  if (p == 1) return(1)
  (p + (p - 2) * rho^2 - 2 * (p - 1) * rho) / (1 - rho^2)
}

# Monte Carlo estimation of theoretical LDA error rate.
# For each grid point, samples from mvrnorm, computes sample pooled-inverse,
# then evaluates lda_error_rate with the true sigma.
simulate_theoretical_error <- function(nk_seq, p_seq, make_sigma, delta, reps, seed) {
  set.seed(seed)
  out_gamma <- numeric(length(nk_seq))
  out_error <- numeric(length(nk_seq))

  for (i in seq_along(nk_seq)) {
    n <- as.integer(nk_seq[[i]])
    p <- max(1L, as.integer(p_seq[[i]]))
    mu1 <- rep(delta, p)
    mu2 <- rep(0, p)
    sigma <- make_sigma(p)

    errs <- numeric(reps)
    for (r in seq_len(reps)) {
      g1 <- MASS::mvrnorm(n, mu1, sigma)
      g2 <- MASS::mvrnorm(n, mu2, sigma)
      xbar1 <- if (p == 1) mean(g1) else colMeans(g1)
      xbar2 <- if (p == 1) mean(g2) else colMeans(g2)
      if (p == 1) {
        g1 <- matrix(g1, ncol = 1)
        g2 <- matrix(g2, ncol = 1)
      }
      s_pooled <- (stats::cov(g1) + stats::cov(g2)) / 2
      s_inv <- tryCatch(ginv_svd(s_pooled), error = function(e) ginv_symmetric(s_pooled))
      errs[[r]] <- lda_error_rate(mu1, mu2, sigma, xbar1, xbar2, s_inv)
    }

    out_gamma[[i]] <- p / (2 * n)
    out_error[[i]] <- mean(errs)
    message(sprintf(
      "theory %02d/%02d: n=%d p=%d gamma=%.3f error=%.6f",
      i, length(nk_seq), n, p, out_gamma[[i]], out_error[[i]]
    ))
  }

  data.frame(gamma = out_gamma, error = out_error)
}

if (should_run("theory")) {
  theory_outputs <- c(
    theoretical_compoundsymmetric = file.path(inputs_dir, "theoretical_compoundsymmetric.csv"),
    theoretical_ar1 = file.path(inputs_dir, "theoretical_ar1.csv"),
    wang_identity_delta = file.path(inputs_dir, "wang_identity_delta.csv"),
    wang_compoundsymmetric = file.path(inputs_dir, "wang_compoundsymmetric.csv"),
    wang_ar1 = file.path(inputs_dir, "wang_ar1.csv")
  )

  # Wang error (analytic)
  wang_nk <- floor(seq(50, 200, length.out = 20))
  wang_gamma_seq <- seq(0.001, 1, length.out = 20)
  wang_p <- ceiling(wang_nk * wang_gamma_seq)
  wang_gamma <- wang_p / wang_nk

  # wang_identity_delta.csv: delta = {0.5, 1, 2}
  wang_id_path <- theory_outputs[["wang_identity_delta"]]
  if (force || !file.exists(wang_id_path)) {
    message("\nComputing Wang identity-delta (analytic)")
    deltas <- c(0.5, 1, 2)
    cols <- lapply(deltas, function(d) {
      maha <- wang_p * d^2
      err <- wang_error_equal_n(maha, wang_gamma)
      data.frame(gamma = wang_gamma, error = err)
    })
    write_csv(do.call(cbind, cols), wang_id_path)
    message("Wrote: ", wang_id_path)
  }

  # wang_compoundsymmetric.csv: delta = 1, rho = {0.1, 0.3, 0.5}
  wang_cs_path <- theory_outputs[["wang_compoundsymmetric"]]
  if (force || !file.exists(wang_cs_path)) {
    message("\nComputing Wang compound-symmetric (analytic)")
    rhos <- c(0.1, 0.3, 0.5)
    cols <- lapply(rhos, function(rho) {
      maha <- wang_p / (1 - rho + rho * wang_p)
      err <- wang_error_equal_n(maha, wang_gamma)
      data.frame(gamma = wang_gamma, error = err)
    })
    write_csv(do.call(cbind, cols), wang_cs_path)
    message("Wrote: ", wang_cs_path)
  }

  # wang_ar1.csv: delta = 1, rho = {0.1, 0.3, 0.5}
  wang_ar1_path <- theory_outputs[["wang_ar1"]]
  if (force || !file.exists(wang_ar1_path)) {
    message("\nComputing Wang AR(1) (analytic)")
    rhos <- c(0.1, 0.3, 0.5)
    cols <- lapply(rhos, function(rho) {
      maha <- vapply(wang_p, function(p) ar1_ones_maha(p, rho), numeric(1))
      err <- wang_error_equal_n(maha, wang_gamma)
      data.frame(gamma = wang_gamma, error = err)
    })
    write_csv(do.call(cbind, cols), wang_ar1_path)
    message("Wrote: ", wang_ar1_path)
  }

  # Theoretical error (Monte Carlo)
  the_nk <- floor(seq(50, 100, length.out = 20))
  the_gamma_seq <- c(seq(0.01, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  the_p <- ceiling(the_nk * 2 * the_gamma_seq)

  the_meta <- function(kind, params, seed_offset) {
    list(
      version = pipeline_version,
      kind = kind,
      params = params,
      n_per_class = the_nk,
      p_dim = the_p,
      reps = reps_theory,
      seed = as.integer(seed + seed_offset)
    )
  }

  # theoretical_compoundsymmetric.csv: rho = {0.1, 0.3, 0.5}
  the_cs_path <- theory_outputs[["theoretical_compoundsymmetric"]]
  the_cs_meta <- the_meta("theoretical_cs", list(delta = 1, rhos = c(0.1, 0.3, 0.5)), 400L)
  if (need_output(the_cs_path, the_cs_meta)) {
    message("\nGenerating theoretical compound-symmetric (Monte Carlo)")
    rhos <- c(0.1, 0.3, 0.5)
    cols <- lapply(seq_along(rhos), function(k) {
      simulate_theoretical_error(
        the_nk, the_p,
        make_sigma = function(p) cov_compound_symmetric(p, rhos[[k]]),
        delta = 1, reps = reps_theory, seed = seed + 400L + k
      )
    })
    write_csv(do.call(cbind, cols), the_cs_path)
    saveRDS(the_cs_meta, meta_path_for(the_cs_path))
    message("Wrote: ", the_cs_path)
  }

  # theoretical_ar1.csv: rho = {0.3, 0.5, 0.7}
  the_ar1_path <- theory_outputs[["theoretical_ar1"]]
  the_ar1_meta <- the_meta("theoretical_ar1", list(delta = 1, rhos = c(0.3, 0.5, 0.7)), 410L)
  if (need_output(the_ar1_path, the_ar1_meta)) {
    message("\nGenerating theoretical AR(1) (Monte Carlo)")
    rhos <- c(0.3, 0.5, 0.7)
    cols <- lapply(seq_along(rhos), function(k) {
      simulate_theoretical_error(
        the_nk, the_p,
        make_sigma = function(p) cov_ar1(p, rhos[[k]]),
        delta = 1, reps = reps_theory, seed = seed + 410L + k
      )
    })
    write_csv(do.call(cbind, cols), the_ar1_path)
    saveRDS(the_ar1_meta, meta_path_for(the_ar1_path))
    message("Wrote: ", the_ar1_path)
  }
}

if (should_run("identity")) {
  identity_out <- file.path(inputs_dir, "identity_empirical_gamma0to5.csv")
  {
  nk_seq <- floor(seq(50, 100, length.out = 20))
  gamma_seq <- c(seq(0.01, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)

	  meta <- list(
	    version = pipeline_version,
	    kind = "gaussian_identity_empirical_delta",
	    deltas = c(0.5, 1, 1.5),
	    n_per_class = nk_seq,
	    p_dim = p_seq,
	    train_frac = 0.7,
	    reps = reps_identity,
	    seed = as.integer(seed + 100L)
	  )

  if (need_output(identity_out, meta)) {
    message("\nGenerating identity empirical inputs (Monte Carlo)")

    mk_sampler <- function(delta) {
      function(n, p) {
        mu1 <- rep(0, p)
        mu2 <- rep(delta, p)
        list(
          g1 = sample_normal_identity(n, mu1),
          g2 = sample_normal_identity(n, mu2)
        )
      }
    }

	    sim05 <- simulate_lda_curve_train_test(nk_seq, p_seq, mk_sampler(0.5), reps = reps_identity, seed = seed + 101L)
	    sim1 <- simulate_lda_curve_train_test(nk_seq, p_seq, mk_sampler(1.0), reps = reps_identity, seed = seed + 102L)
	    sim15 <- simulate_lda_curve_train_test(nk_seq, p_seq, mk_sampler(1.5), reps = reps_identity, seed = seed + 103L)

    # Dummy first row at gamma=1.
    dummy <- function(test_err) data.frame(train_error = 0, test_error = test_err, gamma_result = 1)
    sim05 <- rbind(dummy(0.3608333), sim05)
    sim1 <- rbind(dummy(0.2895833), sim1)
    sim15 <- rbind(dummy(0.114), sim15)

    out_df <- cbind(sim05, sim1, sim15)
    write_csv(out_df, identity_out)
    saveRDS(meta, meta_path_for(identity_out))
    message("Wrote: ", identity_out)
  }
  }
}

if (should_run("cs")) {
cs_out <- file.path(inputs_dir, "compoundsymmetric_empirical_new.csv")
{
  nk_seq <- floor(seq(50, 100, length.out = 20))
  gamma_seq <- c(seq(0.01, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)
  delta <- 1
  rhos <- c(0.1, 0.3, 0.5)

	  meta <- list(
	    version = pipeline_version,
	    kind = "gaussian_compound_symmetric_empirical",
	    delta = delta,
	    rhos = rhos,
	    n_per_class = nk_seq,
	    p_dim = p_seq,
	    train_frac = 0.7,
	    reps = reps_cs,
	    seed = as.integer(seed + 110L)
	  )

  if (need_output(cs_out, meta)) {
    message("\nGenerating compound-symmetric empirical inputs (Monte Carlo)")

    mk_sampler <- function(rho) {
      function(n, p) {
        mu1 <- rep(0, p)
        mu2 <- rep(delta, p)
        list(
          g1 = sample_normal_compound_symmetric(n, mu1, rho),
          g2 = sample_normal_compound_symmetric(n, mu2, rho)
        )
      }
    }

	    sim_r1 <- simulate_lda_curve_train_test(nk_seq, p_seq, mk_sampler(rhos[[1]]), reps = reps_cs, seed = seed + 111L)
	    sim_r3 <- simulate_lda_curve_train_test(nk_seq, p_seq, mk_sampler(rhos[[2]]), reps = reps_cs, seed = seed + 112L)
	    sim_r5 <- simulate_lda_curve_train_test(nk_seq, p_seq, mk_sampler(rhos[[3]]), reps = reps_cs, seed = seed + 113L)

    out_df <- cbind(sim_r1, sim_r3, sim_r5)
    write_csv(out_df, cs_out)
    saveRDS(meta, meta_path_for(cs_out))
    message("Wrote: ", cs_out)
  }
}
}

if (should_run("ar1")) {
ar1_out <- file.path(inputs_dir, "ar1_empirical_new.csv")
{
  nk_seq <- floor(seq(20, 50, length.out = 20))
  gamma_seq <- c(seq(0.01, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)
  delta <- 1
  rhos <- c(0.3, 0.5, 0.7)

	  meta <- list(
	    version = pipeline_version,
	    kind = "gaussian_ar1_empirical",
	    delta = delta,
	    rhos = rhos,
	    n_per_class = nk_seq,
	    p_dim = p_seq,
	    train_frac = 0.7,
	    reps = reps_ar1,
	    seed = as.integer(seed + 120L)
	  )

  if (need_output(ar1_out, meta)) {
    message("\nGenerating AR(1) empirical inputs (Monte Carlo)")

    mk_sampler <- function(rho) {
      function(n, p) {
        mu1 <- rep(0, p)
        mu2 <- rep(delta, p)
        list(
          g1 = sample_normal_ar1(n, mu1, rho),
          g2 = sample_normal_ar1(n, mu2, rho)
        )
      }
    }

	    sim03 <- simulate_lda_curve(nk_seq, p_seq, mk_sampler(rhos[[1]]), reps = reps_ar1, seed = seed + 121L)
	    sim05 <- simulate_lda_curve(nk_seq, p_seq, mk_sampler(rhos[[2]]), reps = reps_ar1, seed = seed + 122L)
	    sim07 <- simulate_lda_curve(nk_seq, p_seq, mk_sampler(rhos[[3]]), reps = reps_ar1, seed = seed + 123L)

    out_df <- data.frame(
      gamma_03 = sim03[["gamma_result"]], err_03 = sim03[["result_sim1"]],
      gamma_05 = sim05[["gamma_result"]], err_05 = sim05[["result_sim1"]],
      gamma_07 = sim07[["gamma_result"]], err_07 = sim07[["result_sim1"]]
    )
    write_csv(out_df, ar1_out)
    saveRDS(meta, meta_path_for(ar1_out))
    message("Wrote: ", ar1_out)
  }
}
}

if (should_run("uniform")) {
unif_files <- file.path(inputs_dir, c("unif_05_025075.csv", "unif_007.csv", "unif_0106.csv"))
{
  nk_seq <- floor(seq(50, 100, length.out = 30))[1:20]
  gamma_seq <- c(seq(0.05, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)

	  meta_unif <- function(min2, max2, seed_offset) {
	    list(
	      version = pipeline_version,
	      kind = "uniform",
	      group2 = c(min = min2, max = max2),
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_unif,
	      seed = as.integer(seed + seed_offset)
	    )
	  }

  need_any_unif <- need_output(unif_files[[1]], meta_unif(0.25, 0.75, 10L)) ||
    need_output(unif_files[[2]], meta_unif(0.00, 0.70, 11L)) ||
    need_output(unif_files[[3]], meta_unif(0.10, 0.60, 12L))

  if (need_any_unif) message("\nGenerating Uniform-distribution inputs (Monte Carlo)")

	  gen_unif <- function(min2, max2, out_path, seed_offset) {
	    meta <- list(
	      version = pipeline_version,
	      kind = "uniform",
	      group2 = c(min = min2, max = max2),
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_unif,
	      seed = as.integer(seed + seed_offset)
	    )
    if (!need_output(out_path, meta)) return(invisible(NULL))
    message("-> ", basename(out_path), sprintf(" (group2 ~ Unif(%.2f, %.2f))", min2, max2))
	    df <- simulate_lda_curve(
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      sampler_pair = function(n, p) list(
	        g1 = sample_uniform(n, p, 0, 0.5),
	        g2 = sample_uniform(n, p, min2, max2)
	      ),
	      reps = reps_unif,
	      seed = seed + seed_offset
	    )
    write_csv(df, out_path)
    saveRDS(meta, meta_path_for(out_path))
    message("Wrote: ", out_path)
  }

  gen_unif(0.25, 0.75, unif_files[[1]], seed_offset = 10L)
  gen_unif(0.00, 0.70, unif_files[[2]], seed_offset = 11L)
  gen_unif(0.10, 0.60, unif_files[[3]], seed_offset = 12L)
}
}

if (should_run("signal")) {
signal_files <- file.path(
  inputs_dir,
  c(
    "signal_1_15.csv", "signal_4_15.csv", "signal_6_15.csv",
    "signal_4_10.csv", "signal_4_20.csv", "signal_4_30.csv"
  )
)

	{
		  delta <- 1
		  nk_seq <- floor(seq(20, 50, length.out = 20))
		  gamma_seq <- c(seq(0.1, 1, length.out = 10), seq(1.02, 4, length.out = 10))
		  p_seq <- ceiling(nk_seq * 2 * gamma_seq)

		  reps_noise <- resolve_reps(100L)
		  reps_redund_10_20 <- resolve_reps(80L)
		  reps_redund_30 <- resolve_reps(100L)

		  # Additive Gaussian noise on a fraction of columns.
				  meta_noise <- function(noise_sd, noise_frac, seed_offset, reps) {
				    list(
				      version = pipeline_version,
				      kind = "signal_noise",
				      noise_sd = noise_sd,
				      noise_frac = noise_frac,
				      noise_col_rounding = "ceiling_(p+1)",
				      noise_apply = "additive_gaussian",
				      seed_mode = "offset",
				      delta = delta,
				      base_sd = 1,
			      n_per_class = nk_seq,
		      p_dim = p_seq,
		      train_frac = 0.7,
		      reps = as.integer(reps),
		      seed = as.integer(seed + seed_offset)
		    )
		  }

		  # Redundant features: a fraction of columns zeroed out.
				  meta_redund <- function(redundant_frac, seed_offset, reps) {
				    list(
				      version = pipeline_version,
				      kind = "redundant_features",
				      redundant_frac = redundant_frac,
				      redundant_mode = "zero_columns",
				      redundant_col_rounding = "ceiling_(p+1)",
				      seed_mode = "offset",
				      delta = delta,
				      base_sd = 1,
			      n_per_class = nk_seq,
		      p_dim = p_seq,
		      train_frac = 0.7,
		      reps = as.integer(reps),
		      seed = as.integer(seed + seed_offset)
		    )
		  }

		  need_any_signal <- need_output(file.path(inputs_dir, "signal_1_15.csv"), meta_noise(1, 0.15, 0L, reps_noise)) ||
		    need_output(file.path(inputs_dir, "signal_4_15.csv"), meta_noise(4, 0.15, 0L, reps_noise)) ||
		    need_output(file.path(inputs_dir, "signal_6_15.csv"), meta_noise(6, 0.15, 0L, reps_noise)) ||
		    need_output(file.path(inputs_dir, "signal_4_10.csv"), meta_redund(0.10, 0L, reps_redund_10_20)) ||
		    need_output(file.path(inputs_dir, "signal_4_20.csv"), meta_redund(0.20, 0L, reps_redund_10_20)) ||
		    need_output(file.path(inputs_dir, "signal_4_30.csv"), meta_redund(0.30, 0L, reps_redund_30))

  if (need_any_signal) message("\nGenerating signal/noise + redundant-feature inputs (Monte Carlo)")

		  gen_signal_noise <- function(noise_sd, noise_frac, out_path, seed_offset, reps) {
		    meta <- meta_noise(noise_sd = noise_sd, noise_frac = noise_frac, seed_offset = seed_offset, reps = reps)
		    if (!need_output(out_path, meta)) return(invisible(NULL))
		    message("-> ", basename(out_path), sprintf(" (noise sd=%.1f, noisy cols=%.0f%%)", noise_sd, 100 * noise_frac))

					    df <- simulate_lda_curve(
					      n_per_class = nk_seq,
					      p_dim = p_seq,
				      sampler_pair = function(n, p) {
	        mu1 <- numeric(p)
	        mu2 <- rep(delta, p)
	        g1 <- matrix(stats::rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p) + rep(mu1, each = n)
	        g2 <- matrix(stats::rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p) + rep(mu2, each = n)

			        n_noisy <- max(1L, as.integer(ceiling(noise_frac * (p + 1L))))
			        n_noisy <- min(p, n_noisy)
		        noisy_cols <- sample.int(p, size = n_noisy, replace = FALSE)
			        if (n_noisy > 0) {
		          eps <- matrix(stats::rnorm(2L * n * n_noisy, mean = 0, sd = noise_sd), nrow = 2L * n, ncol = n_noisy)
		          g1[, noisy_cols] <- g1[, noisy_cols, drop = FALSE] + eps[seq_len(n), , drop = FALSE]
		          g2[, noisy_cols] <- g2[, noisy_cols, drop = FALSE] + eps[seq.int(n + 1L, 2L * n), , drop = FALSE]
		        }

		        list(g1 = g1, g2 = g2)
						      },
						      reps = reps,
						      seed = seed + seed_offset,
						      fit_fn = lda_ginv
						    )
				    write_csv(df, out_path)
				    saveRDS(meta, meta_path_for(out_path))
				    message("Wrote: ", out_path)
				  }

		  gen_redundant_features <- function(redundant_frac, out_path, seed_offset, reps) {
		    meta <- meta_redund(redundant_frac = redundant_frac, seed_offset = seed_offset, reps = reps)
		    if (!need_output(out_path, meta)) return(invisible(NULL))
		    message("-> ", basename(out_path), sprintf(" (irrelevant cols=%.0f%%)", 100 * redundant_frac))

					    df <- simulate_lda_curve(
					      n_per_class = nk_seq,
					      p_dim = p_seq,
				      sampler_pair = function(n, p) {
		        mu1 <- numeric(p)
		        mu2 <- rep(delta, p)
			        g1 <- matrix(stats::rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p) + rep(mu1, each = n)
			        g2 <- matrix(stats::rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p) + rep(mu2, each = n)

			        n_redundant <- max(1L, as.integer(ceiling(redundant_frac * (p + 1L))))
			        n_redundant <- min(p, n_redundant)
			        redundant_cols <- sample.int(p, size = n_redundant, replace = FALSE)
		        if (n_redundant > 0) {
		          g1[, redundant_cols] <- 0
	          g2[, redundant_cols] <- 0
	        }
			        list(g1 = g1, g2 = g2)
						      },
						      reps = reps,
						      seed = seed + seed_offset,
						      fit_fn = lda_ginv
						    )
				    write_csv(df, out_path)
				    saveRDS(meta, meta_path_for(out_path))
				    message("Wrote: ", out_path)
				  }

	  # Varying noise level.
		  gen_signal_noise(1, 0.15, file.path(inputs_dir, "signal_1_15.csv"), seed_offset = 0L, reps = reps_noise)
		  gen_signal_noise(4, 0.15, file.path(inputs_dir, "signal_4_15.csv"), seed_offset = 0L, reps = reps_noise)
		  gen_signal_noise(6, 0.15, file.path(inputs_dir, "signal_6_15.csv"), seed_offset = 0L, reps = reps_noise)

	  # Varying redundant feature percentage.
		  gen_redundant_features(0.10, file.path(inputs_dir, "signal_4_10.csv"), seed_offset = 0L, reps = reps_redund_10_20)
		  gen_redundant_features(0.20, file.path(inputs_dir, "signal_4_20.csv"), seed_offset = 0L, reps = reps_redund_10_20)
		  gen_redundant_features(0.30, file.path(inputs_dir, "signal_4_30.csv"), seed_offset = 0L, reps = reps_redund_30)
}
}

if (should_run("tdist_id")) {
tdist_files <- file.path(inputs_dir, c("tdist_id_df1.csv", "tdist_id_df4.csv", "tdist_id_df30.csv"))
{
  nk_seq <- floor(seq(50, 100, length.out = 20))
  gamma_seq <- c(seq(0.05, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)

	  meta_tdist <- function(df, seed_offset) {
	    list(
	      version = pipeline_version,
	      kind = "t_identity",
	      df = df,
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_tdist_id,
	      seed = as.integer(seed + seed_offset)
	    )
	  }

  need_any_tdist <- need_output(tdist_files[[1]], meta_tdist(1, 30L)) ||
    need_output(tdist_files[[2]], meta_tdist(4, 31L)) ||
    need_output(tdist_files[[3]], meta_tdist(30, 32L))

  if (need_any_tdist) message("\nGenerating t-distribution (identity) inputs (Monte Carlo)")

  sample_t_id <- function(n, p, mu, df) {
    z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    w <- stats::rchisq(n, df = df) / df
    sweep(z, 1, sqrt(w), "/") + rep(mu, each = n)
  }

	  gen_tdist <- function(df, out_path, seed_offset) {
	    meta <- list(
	      version = pipeline_version,
	      kind = "t_identity",
	      df = df,
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_tdist_id,
	      seed = as.integer(seed + seed_offset)
	    )
    if (!need_output(out_path, meta)) return(invisible(NULL))
    message("-> ", basename(out_path), sprintf(" (df=%d)", df))

	    df_out <- simulate_lda_curve(
	      n_per_class = nk_seq,
	      p_dim = p_seq,
      sampler_pair = function(n, p) {
        mu1 <- rep(0, p)
        mu2 <- rep(1, p)
        list(
          g1 = sample_t_id(n, p, mu = mu1, df = df),
          g2 = sample_t_id(n, p, mu = mu2, df = df)
        )
      },
	      reps = reps_tdist_id,
	      seed = seed + seed_offset
	    )
    write_csv(df_out, out_path)
    saveRDS(meta, meta_path_for(out_path))
    message("Wrote: ", out_path)
  }

  gen_tdist(1, tdist_files[[1]], seed_offset = 30L)
  gen_tdist(4, tdist_files[[2]], seed_offset = 31L)
  gen_tdist(30, tdist_files[[3]], seed_offset = 32L)
}
}

if (should_run("poisson")) {
pois_specs <- list(
  poisson_unif_23_34 = c(2, 3, 3, 4),
  poisson_unif_24_46 = c(2, 4, 4, 6),
  poisson_unif_67_78 = c(6, 7, 7, 8)
)

{
  nk_seq <- floor(seq(50, 100, length.out = 30))
  gamma_seq_base <- c(seq(0.01, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  gamma_seq <- rep(gamma_seq_base, length.out = length(nk_seq))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)

  for (nm in names(pois_specs)) {
    out_path <- file.path(inputs_dir, paste0(nm, ".csv"))
    spec <- pois_specs[[nm]]

	    meta <- list(
	      version = pipeline_version,
	      kind = "poisson",
      lambda1 = c(min = spec[[1]], max = spec[[2]]),
      lambda2 = c(min = spec[[3]], max = spec[[4]]),
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_poisson,
	      seed = as.integer(seed + 200L)
	    )

    if (!need_output(out_path, meta)) next
    message("\nGenerating ", nm, " (Monte Carlo)")

	    df <- simulate_lda_curve(
	      n_per_class = nk_seq,
	      p_dim = p_seq,
      sampler_pair = function(n, p) list(
        g1 = sample_multi_poisson(n, p, lambda_min = spec[[1]], lambda_max = spec[[2]]),
        g2 = sample_multi_poisson(n, p, lambda_min = spec[[3]], lambda_max = spec[[4]])
      ),
	      reps = reps_poisson,
	      seed = seed + 200L
	    )
    write_csv(df, out_path)
    saveRDS(meta, meta_path_for(out_path))
    message("Wrote: ", out_path)
  }
}
}

if (should_run("tdist_cs")) {
cs_tdist_files <- file.path(inputs_dir, c("cs_df1_03.csv", "cs_df4_03.csv", "cs_df30_03.csv"))
{
  nk_seq <- floor(seq(50, 100, length.out = 20))
  gamma_seq <- c(seq(0.05, 1, length.out = 10), seq(1.02, 3, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)
  rho <- 0.3

	  meta_cs_t <- function(df, seed_offset) {
	    list(
	      version = pipeline_version,
	      kind = "t_compound_symmetric",
      rho = rho,
      df = df,
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_tdist_cs,
	      seed = as.integer(seed + seed_offset)
	    )
	  }

  gen_cs_t <- function(df, out_path, seed_offset) {
    meta <- meta_cs_t(df, seed_offset)
    if (!need_output(out_path, meta)) return(invisible(NULL))
    message("\nGenerating ", basename(out_path), sprintf(" (t, CS rho=%.1f, df=%d)", rho, df))

	    df_out <- simulate_lda_curve(
	      n_per_class = nk_seq,
	      p_dim = p_seq,
      sampler_pair = function(n, p) {
        mu1 <- rep(0, p)
        mu2 <- rep(1, p)
        list(
          g1 = rmvt_compound_symmetric(n, mu = mu1, rho = rho, df = df),
          g2 = rmvt_compound_symmetric(n, mu = mu2, rho = rho, df = df)
        )
      },
	      reps = reps_tdist_cs,
	      seed = seed + seed_offset
	    )
    write_csv(df_out, out_path)
    saveRDS(meta, meta_path_for(out_path))
    message("Wrote: ", out_path)
  }

  gen_cs_t(1, cs_tdist_files[[1]], seed_offset = 210L)
  gen_cs_t(4, cs_tdist_files[[2]], seed_offset = 211L)
  gen_cs_t(30, cs_tdist_files[[3]], seed_offset = 212L)
}
}

if (should_run("flip")) {
flip_specs <- list(
  flip_05 = 0.05,
  flip_10 = 0.10,
  flip_20 = 0.20
)

{
  nk_seq <- floor(seq(20, 50, length.out = 20))
  gamma_seq <- c(seq(0.1, 1, length.out = 10), seq(1.02, 4, length.out = 10))
  p_seq <- ceiling(nk_seq * 2 * gamma_seq)
  delta <- 1

  for (nm in names(flip_specs)) {
    out_path <- file.path(inputs_dir, paste0(nm, ".csv"))
    flip_frac <- flip_specs[[nm]]

	    meta <- list(
	      version = pipeline_version,
	      kind = "label_flip",
      delta = delta,
      flip_frac = flip_frac,
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      train_frac = 0.7,
	      reps = reps_flip,
	      seed = as.integer(seed + 300L)
	    )

    if (!need_output(out_path, meta)) next
    message("\nGenerating ", nm, " (Monte Carlo)")

	    df <- simulate_lda_curve_label_flip(
	      n_per_class = nk_seq,
	      p_dim = p_seq,
	      delta = delta,
	      flip_frac = flip_frac,
	      reps = reps_flip,
	      seed = seed + 300L
	    )
    write_csv(df, out_path)
    saveRDS(meta, meta_path_for(out_path))
    message("Wrote: ", out_path)
  }
}
}

message("\nDone. Inputs live in: ", inputs_dir)
