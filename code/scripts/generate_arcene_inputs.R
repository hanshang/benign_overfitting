#!/usr/bin/env Rscript

# Regenerate ARCENE-related CSV inputs in code/paper/inputs/.

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
arcene_dir <- file.path(code_root, "ARCENE")

dir.create(inputs_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(script_dir, "experiments_lib.R"))

args <- commandArgs(trailingOnly = TRUE)
has_flag <- function(flag) any(args == flag)
get_arg_value <- function(prefix) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit) == 0) return(NULL)
  sub(paste0("^", prefix), "", hit[[1]])
}

force <- has_flag("--force")

# Parse --seed (default 123).
seed <- 123L
seed_arg <- get_arg_value("--seed=")
if (!is.null(seed_arg)) {
  seed <- suppressWarnings(as.integer(seed_arg))
  if (is.na(seed)) stop("Invalid --seed value: ", seed_arg)
}

# Parse --trials / --arcene-trials (default 30).
trials <- 30L
trials_arg <- get_arg_value("--trials=")
if (is.null(trials_arg)) trials_arg <- get_arg_value("--arcene-trials=")
if (!is.null(trials_arg)) {
  trials <- suppressWarnings(as.integer(trials_arg))
  if (is.na(trials) || trials < 1) stop("Invalid --trials value: ", trials_arg)
}

# Stop if file does not exist.
assert_exists <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
  invisible(path)
}

# Read one ARCENE split (data + labels) into list(x, y).
read_arcene_split <- function(data_path, label_path, drop_first_row = TRUE) {
  assert_exists(data_path, "ARCENE data file")
  assert_exists(label_path, "ARCENE label file")

  y_full <- scan(label_path, quiet = TRUE)
  x_raw <- scan(data_path, quiet = TRUE, what = double())
  if (length(x_raw) %% length(y_full) != 0) {
    stop("ARCENE data length not divisible by label count: ", data_path)
  }
  x_full <- matrix(x_raw, nrow = length(y_full), byrow = TRUE)

  if (!drop_first_row) return(list(x = x_full, y = y_full))

  list(x = x_full[-1, , drop = FALSE], y = y_full[-1])
}

# Load ARCENE train + validation splits.
load_arcene <- function(drop_first_row = TRUE) {
  tr <- read_arcene_split(
    file.path(arcene_dir, "arcene_train.data"),
    file.path(arcene_dir, "arcene_train.labels"),
    drop_first_row = drop_first_row
  )
  te <- read_arcene_split(
    file.path(arcene_dir, "arcene_valid.data"),
    file.path(arcene_dir, "arcene_valid.labels"),
    drop_first_row = drop_first_row
  )

  if (nrow(tr$x) != length(tr$y)) stop("ARCENE train x/y mismatch.")
  if (nrow(te$x) != length(te$y)) stop("ARCENE valid x/y mismatch.")

  list(train = tr, test = te)
}

# Drop first column from train and test feature matrices.
preprocess_arcene_features <- function(x_train, x_test) {
  if (ncol(x_train) < 2) stop("ARCENE preprocessing removed too many features.")

  # Drop the first column (paper convention).
  x_train <- x_train[, -1, drop = FALSE]
  x_test <- x_test[, -1, drop = FALSE]
  list(train = x_train, test = x_test)
}

# Predict {-1, 1} from a linear classifier (beta, thresh).
predict_linear <- function(x, beta, thresh) {
  score <- as.vector(x %*% beta)
  ifelse(score > as.numeric(thresh), 1, -1)
}

# Diagonal LDA with pooled per-feature variances.
dlda_fit_pooled <- function(x, y, tiny_constant = 1e-2) {
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible sizes.")

  x1 <- x[y == 1, , drop = FALSE]
  x2 <- x[y == -1, , drop = FALSE]
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  if (n1 < 2 || n2 < 2) stop("Need at least 2 samples per class.")

  mu1 <- colMeans(x1)
  mu2 <- colMeans(x2)

  v1 <- apply(x1, 2, stats::var)
  v2 <- apply(x2, 2, stats::var)
  v <- ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
  v[v == 0] <- tiny_constant

  beta <- (mu1 - mu2) / v
  thresh <- 0.5 * sum(beta * (mu1 + mu2)) - log(n1 / n2)
  list(beta = beta, thresh = thresh, mu1 = mu1, mu2 = mu2, v = v, n1 = n1, n2 = n2)
}

# Shrinkage LDA via Ledoit-Wolf correlation shrinkage.
slda_fit_shrunk_corr <- function(x, y, tiny_constant = 1e-2) {
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible sizes.")

  x1 <- x[y == 1, , drop = FALSE]
  x2 <- x[y == -1, , drop = FALSE]
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  n <- n1 + n2
  df <- n - 2
  if (n1 < 2 || n2 < 2) stop("Need at least 2 samples per class.")

  mu1 <- colMeans(x1)
  mu2 <- colMeans(x2)
  e <- rbind(sweep(x1, 2, mu1, "-"), sweep(x2, 2, mu2, "-"))

  v <- colSums(e^2) / df
  v[v == 0] <- tiny_constant

  u <- sweep(e, 2, sqrt(v), "/")

  denom <- n
  r <- crossprod(u) / denom
  rho <- sum(r^2) - sum(diag(r)^2)

  row_sums2 <- rowSums(u^2)
  row_sums4 <- rowSums(u^4)
  a_sum <- sum(row_sums2^2 - row_sums4)

  phi_off <- (a_sum / denom - rho) / denom
  lambda <- if (rho > 0) max(0, min(1, phi_off / rho)) else 1

  dif_c <- (mu1 - mu2) / sqrt(v)

  # Solve R^-1 dif_c via SVD.
  z <- u / sqrt(denom)
  sv <- base::svd(z, nu = 0)
  d <- sv$d
  vmat <- sv$v
  eig <- d^2

  if (lambda == 0) {
    inv_eig <- ifelse(eig > 0, 1 / eig, 0)
    beta_c <- vmat %*% (inv_eig * as.vector(crossprod(vmat, dif_c)))
  } else {
    base_inv <- 1 / lambda
    inv_shr_eig <- 1 / ((1 - lambda) * eig + lambda)
    proj <- as.vector(crossprod(vmat, dif_c))
    beta_c <- base_inv * dif_c + vmat %*% ((inv_shr_eig - base_inv) * proj)
  }

  beta <- as.vector(beta_c / sqrt(v))
  thresh <- 0.5 * sum(beta * (mu1 + mu2)) - log(n1 / n2)
  list(beta = beta, thresh = thresh, lambda = lambda)
}

# Modified LDA with eigenvalue flooring on the correlation matrix.
mlda_fit_floor_corr <- function(x, y, floor_eig = 1, tiny_constant = 1e-2) {
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible sizes.")

  x1 <- x[y == 1, , drop = FALSE]
  x2 <- x[y == -1, , drop = FALSE]
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  n <- n1 + n2
  df <- n - 2
  if (n1 < 2 || n2 < 2) stop("Need at least 2 samples per class.")

  mu1 <- colMeans(x1)
  mu2 <- colMeans(x2)
  e <- rbind(sweep(x1, 2, mu1, "-"), sweep(x2, 2, mu2, "-"))

  w0 <- e / sqrt(df)
  v <- colSums(w0^2)
  v[v == 0] <- tiny_constant

  z <- sweep(w0, 2, sqrt(v), "/")

  sv <- base::svd(z, nu = 0)
  d <- sv$d
  vmat <- sv$v
  eig <- d^2

  inv_eig <- 1 / pmax(eig, floor_eig)

  dif_c <- (mu1 - mu2) / sqrt(v)
  proj <- as.vector(crossprod(vmat, dif_c))
  beta_c <- (1 / floor_eig) * dif_c + vmat %*% ((inv_eig - 1 / floor_eig) * proj)

  beta <- as.vector(beta_c / sqrt(v))
  thresh <- 0.5 * sum(beta * (mu1 + mu2)) - log(n1 / n2)
  list(beta = beta, thresh = thresh)
}

# Feature-count grid for random feature sampling experiment.
arcene_random_p_seq <- function() {
  c(seq(2, 97, by = 5), floor(seq(101, 1000, length.out = 100)))
}

# LDA error vs gamma using random feature subsets.
generate_arcene_random <- function(train, test, seed, trials) {
  x_tr <- train$x
  y_tr <- train$y
  x_te <- test$x
  y_te <- test$y

  n_test <- length(y_te)
  p_seq <- arcene_random_p_seq()
  max_p <- max(p_seq)
  gamma <- p_seq / 99

  total_misclass <- integer(length(p_seq))

  set.seed(seed)
  for (s in seq_len(trials)) {
    idx_full <- sample.int(ncol(x_tr), size = max_p, replace = FALSE)
    for (i in seq_along(p_seq)) {
      p <- p_seq[[i]]
      idx <- idx_full[seq_len(p)]
      mod <- lda_ginv_svd(x_tr[, idx, drop = FALSE], y_tr)
      pred <- predict_lda(x_te[, idx, drop = FALSE], mod)
      total_misclass[[i]] <- total_misclass[[i]] + sum(pred != y_te)
    }
    if (s %% 5 == 0 || s == trials) message("arcene_random: trial ", s, "/", trials)
  }

  data.frame(gamma = gamma, averages = total_misclass / (trials * n_test))
}

# LDA error vs number of principal components.
generate_arcene_pca <- function(train, test) {
  x_tr <- train$x
  y_tr <- train$y
  x_te <- test$x
  y_te <- test$y

  pca <- stats::prcomp(x_tr, center = TRUE, scale. = FALSE)

  test_scaled <- scale(x_te, center = pca$center, scale = pca$scale)
  test_scores <- test_scaled %*% pca$rotation

  pcs <- 2:98
  errs <- numeric(length(pcs))
  for (k in seq_along(pcs)) {
    i <- pcs[[k]]
    mod <- lda_ginv(pca$x[, 1:i, drop = FALSE], y_tr)
    pred <- predict_lda(test_scores[, 1:i, drop = FALSE], mod)
    errs[[k]] <- mean(pred != y_te)
  }

  out <- data.frame(gamma = pcs / 99, test_error = errs)
  names(out)[[1]] <- ""
  out
}

# Compare Pseudo-LDA, DLDA, SLDA, and MLDA on ARCENE.
generate_arcene_combined_methods <- function(train, test) {
  d <- preprocess_arcene_features(train$x, test$x)
  x_tr_full <- d$train
  x_te_full <- d$test
  y_tr <- train$y
  y_te <- test$y

  p_seq <- floor(seq(2, 1000, length.out = 30))
  gamma <- p_seq / 99

  out_error <- numeric(length(p_seq))
  out_dlda <- numeric(length(p_seq))
  out_slda <- numeric(length(p_seq))
  out_mlda <- numeric(length(p_seq))

  for (i in seq_along(p_seq)) {
    p <- p_seq[[i]]
    x_tr <- x_tr_full[, 1:p, drop = FALSE]
    x_te <- x_te_full[, 1:p, drop = FALSE]

    mod_plda <- lda_ginv(x_tr, y_tr)
    pred_plda <- predict_lda(x_te, mod_plda)
    out_error[[i]] <- mean(pred_plda != y_te)

    mod_dlda <- dlda_fit_pooled(x_tr, y_tr)
    pred_dlda <- predict_linear(x_te, mod_dlda$beta, mod_dlda$thresh)
    out_dlda[[i]] <- mean(pred_dlda != y_te)

    mod_slda <- slda_fit_shrunk_corr(x_tr, y_tr)
    pred_slda <- predict_linear(x_te, mod_slda$beta, mod_slda$thresh)
    out_slda[[i]] <- mean(pred_slda != y_te)

    mod_mlda <- mlda_fit_floor_corr(x_tr, y_tr, floor_eig = 1)
    pred_mlda <- predict_linear(x_te, mod_mlda$beta, mod_mlda$thresh)
    out_mlda[[i]] <- mean(pred_mlda != y_te)

    message(sprintf("ARCENE combined: %02d/%02d p=%d gamma=%.3f", i, length(p_seq), p, gamma[[i]]))
  }

  data.frame(V1 = gamma, error = out_error, dlda = out_dlda, slda = out_slda, mlda = out_mlda)
}

outputs <- list(
  random = file.path(inputs_dir, "arcene_random.csv"),
  pca = file.path(inputs_dir, "arcene_PCA.csv"),
  combined = file.path(inputs_dir, "ARCENE_combined_methods.csv")
)

need_random <- force || !file.exists(outputs$random)
need_pca <- force || !file.exists(outputs$pca)
need_combined <- force || !file.exists(outputs$combined)

if (!need_random && !need_pca && !need_combined) {
  message("ARCENE inputs already present; nothing to do.")
  quit(save = "no", status = 0, runLast = FALSE)
}

message("Loading ARCENE data")
arc <- load_arcene(drop_first_row = TRUE)

if (need_random) {
  message("\nGenerating arcene_random.csv (seed=", seed, ", trials=", trials, ")")
  df <- generate_arcene_random(arc$train, arc$test, seed = seed, trials = trials)
  old_digits <- getOption("digits")
  options(digits = 9)
  on.exit(options(digits = old_digits), add = TRUE)
  utils::write.csv(df, outputs$random, row.names = FALSE, quote = FALSE)
  message("Wrote: ", outputs$random)
}

if (need_combined) {
  message("\nGenerating ARCENE_combined_methods.csv")
  df <- generate_arcene_combined_methods(arc$train, arc$test)
  old_digits <- getOption("digits")
  options(digits = 16)
  on.exit(options(digits = old_digits), add = TRUE)
  utils::write.csv(df, outputs$combined)
  message("Wrote: ", outputs$combined)
}

if (need_pca) {
  message("\nGenerating arcene_PCA.csv")
  df <- generate_arcene_pca(arc$train, arc$test)
  old_digits <- getOption("digits")
  options(digits = 16)
  on.exit(options(digits = old_digits), add = TRUE)
  utils::write.csv(df, outputs$pca, row.names = FALSE)
  message("Wrote: ", outputs$pca)
}
