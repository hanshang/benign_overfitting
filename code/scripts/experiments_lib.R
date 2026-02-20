# Shared functions for LDA simulations.
# Includes LDA fitting/prediction, sampling, covariance constructors, theory formulas.

suppressPackageStartupMessages({
  library(MASS)
})

options(stringsAsFactors = FALSE)

# Stratified train/test split. Returns list(train, test) of indices.
split_train_test <- function(y, train_frac = 0.7, seed = NULL) {
  if (!is.numeric(train_frac) || length(train_frac) != 1 || train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be a scalar in (0, 1).")
  }
  if (!is.null(seed)) set.seed(seed)

  y <- as.vector(y)
  classes <- sort(unique(y))
  idx_train <- integer(0)
  for (cls in classes) {
    idx <- which(y == cls)
    n_train <- ceiling(length(idx) * train_frac)
    n_train <- max(1L, min(length(idx) - 1L, n_train))
    idx_train <- c(idx_train, sample(idx, n_train))
  }
  idx_train <- sort(idx_train)
  idx_test <- setdiff(seq_along(y), idx_train)
  list(train = idx_train, test = idx_test)
}

# Generalized inverse for symmetric matrices.
ginv_symmetric <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x))) {
    stop("'x' must be a numeric or complex matrix")
  }
  if (!is.matrix(x)) x <- as.matrix(x)
  eig <- eigen(x, symmetric = TRUE)
  d <- eig$values
  v <- eig$vectors

  if (length(d) == 0) return(array(0, dim(x)[2L:1L]))
  positive <- d > max(tol * d[1L], 0)

  if (all(positive)) {
    return(v %*% (1 / d * t(v)))
  }
  if (!any(positive)) {
    return(array(0, dim(x)[2L:1L]))
  }
  vpos <- v[, positive, drop = FALSE]
  vpos %*% ((1 / d[positive]) * t(vpos))
}

# Generalized inverse via SVD.
ginv_svd <- function(x, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(x)) > 2L || !(is.numeric(x) || is.complex(x))) {
    stop("'x' must be a numeric or complex matrix")
  }
  if (!is.matrix(x)) x <- as.matrix(x)

  s <- base::svd(x, nu = nrow(x), nv = ncol(x), LINPACK = TRUE)
  d <- s$d
  if (length(d) == 0) return(array(0, dim(x)[2L:1L]))

  positive <- d > max(tol * d[1L], 0)
  if (!any(positive)) return(array(0, dim(x)[2L:1L]))

  u <- s$u[, positive, drop = FALSE]
  v <- s$v[, positive, drop = FALSE]
  v %*% ((1 / d[positive]) * t(u))
}

# LDA via Moore-Penrose pseudoinverse of pooled covariance.
# Returns list(coef, thresh, prec_matrix, group_mean).
lda_ginv <- function(x, y) {
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible sizes.")

  p <- ncol(x)
  mu <- matrix(0, nrow = p, ncol = 2)

  if (p == 1) {
    mu[, 1] <- mean(x[y == 1, , drop = TRUE])
    mu[, 2] <- mean(x[y == -1, , drop = TRUE])
    s <- stats::cov(x)
  } else {
    mu[, 1] <- colMeans(x[y == 1, , drop = FALSE])
    mu[, 2] <- colMeans(x[y == -1, , drop = FALSE])
    s <- (stats::cov(x[y == 1, , drop = FALSE]) + stats::cov(x[y == -1, , drop = FALSE])) / 2
  }

  s_inv <- tryCatch(
    ginv_svd(s),
    error = function(e) ginv_symmetric(s)
  )
  dif_vec <- matrix(mu[, 1] - mu[, 2], ncol = 1)
  sum_vec <- matrix(mu[, 1] + mu[, 2], ncol = 1)

  coef <- s_inv %*% dif_vec
  thresh <- 0.5 * t(coef) %*% sum_vec

  list(coef = coef, thresh = thresh, prec_matrix = s_inv, group_mean = mu)
}

# SVD-based LDA pseudoinverse.
lda_ginv_svd <- function(x, y, tol = NULL) {
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible sizes.")

  p <- ncol(x)
  if (p == 1) return(lda_ginv(x, y))

  x1 <- x[y == 1, , drop = FALSE]
  x2 <- x[y == -1, , drop = FALSE]
  n1 <- nrow(x1)
  n2 <- nrow(x2)
  if (n1 < 2 || n2 < 2) stop("Need at least 2 samples per class for covariance estimation.")

  mu1 <- colMeans(x1)
  mu2 <- colMeans(x2)
  mu <- cbind(mu1, mu2)

  x1c <- sweep(x1, 2, mu1, "-") / sqrt(n1 - 1)
  x2c <- sweep(x2, 2, mu2, "-") / sqrt(n2 - 1)
  w <- rbind(x1c, x2c)

  sv <- base::svd(w, nu = 0)
  d <- sv$d
  v <- sv$v

  if (length(d) == 0) {
    coef <- matrix(0, nrow = p, ncol = 1)
    thresh <- matrix(0, nrow = 1, ncol = 1)
    return(list(coef = coef, thresh = thresh, prec_matrix = NULL, group_mean = mu))
  }

  if (is.null(tol)) {
    tol <- max(dim(w)) * max(d) * .Machine$double.eps
  }
  if (!is.numeric(tol) || length(tol) != 1 || tol < 0) stop("tol must be a non-negative scalar.")

  keep <- d > tol
  if (!any(keep)) {
    coef <- matrix(0, nrow = p, ncol = 1)
    thresh <- matrix(0, nrow = 1, ncol = 1)
    return(list(coef = coef, thresh = thresh, prec_matrix = NULL, group_mean = mu))
  }

  d <- d[keep]
  v <- v[, keep, drop = FALSE]

  dif_vec <- matrix(mu1 - mu2, ncol = 1)
  tmp <- as.vector(crossprod(v, dif_vec)) / (d^2)
  coef <- 2 * v %*% matrix(tmp, ncol = 1)

  sum_vec <- matrix(mu1 + mu2, ncol = 1)
  thresh <- 0.5 * t(coef) %*% sum_vec

  list(coef = coef, thresh = thresh, prec_matrix = NULL, group_mean = mu)
}

# Predict {-1, 1} from an lda_ginv model.
predict_lda <- function(newdata, model) {
  newdata <- as.matrix(newdata)
  score <- as.vector(t(model$coef) %*% t(newdata))
  ifelse(score > as.numeric(model$thresh), 1, -1)
}

# n x p matrix of iid Uniform(min, max) entries.
sample_uniform <- function(n, p, min, max) {
  matrix(stats::runif(n * p, min = min, max = max), nrow = n, ncol = p)
}

# Sample n rows from N(mu, I_p).
sample_normal_identity <- function(n, mu) {
  mu <- as.numeric(mu)
  p <- length(mu)
  matrix(stats::rnorm(n * p), nrow = n, ncol = p) + rep(mu, each = n)
}

# Sample n rows from N(mu, Sigma_CS) where Sigma_ij = rho for i != j.
sample_normal_compound_symmetric <- function(n, mu, rho) {
  if (!is.numeric(rho) || length(rho) != 1 || is.na(rho) || rho < 0 || rho >= 1) {
    stop("rho must be a scalar in [0, 1).")
  }
  mu <- as.numeric(mu)
  p <- length(mu)
  z <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  shared <- stats::rnorm(n)
  x <- sqrt(1 - rho) * z + sqrt(rho) * matrix(shared, nrow = n, ncol = p)
  x + rep(mu, each = n)
}

# Sample n rows from N(mu, Sigma_AR1) where Sigma_ij = rho^|i-j|.
sample_normal_ar1 <- function(n, mu, rho) {
  if (!is.numeric(rho) || length(rho) != 1 || is.na(rho) || abs(rho) >= 1) {
    stop("rho must be a scalar in (-1, 1).")
  }
  mu <- as.numeric(mu)
  p <- length(mu)
  eps <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  x <- matrix(0, nrow = n, ncol = p)
  x[, 1] <- eps[, 1]
  if (p > 1) {
    sd_innov <- sqrt(1 - rho^2)
    for (j in 2:p) x[, j] <- rho * x[, j - 1] + sd_innov * eps[, j]
  }
  x + rep(mu, each = n)
}

# n x p Poisson matrix with per-feature rates lam_j ~ Unif(lambda_min, lambda_max).
sample_multi_poisson <- function(n, p, lambda_min, lambda_max) {
  lambda <- stats::runif(p, min = lambda_min, max = lambda_max)
  out <- matrix(0, nrow = n, ncol = p)
  for (j in seq_len(p)) out[, j] <- stats::rpois(n, lambda[j])
  out
}

# Multivariate t sampler via scale-mixture of normals.
rmvt_base <- function(n, mu, sigma, df) {
  mu <- as.numeric(mu)
  p <- length(mu)
  sigma <- as.matrix(sigma)
  if (!all(dim(sigma) == c(p, p))) stop("sigma must be p x p with p = length(mu).")
  if (!is.numeric(df) || length(df) != 1 || df <= 0) stop("df must be a positive scalar.")

  z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)
  w <- stats::rchisq(n, df = df) / df
  sweep(z, 1, sqrt(w), "/") + rep(mu, each = n)
}

# Multivariate t with compound-symmetric covariance (unit diagonal).
rmvt_compound_symmetric <- function(n, mu, rho, df) {
  mu <- as.numeric(mu)
  p <- length(mu)
  z <- sample_normal_compound_symmetric(n, rep(0, p), rho)
  w <- stats::rchisq(n, df = df) / df
  sweep(z, 1, sqrt(w), "/") + rep(mu, each = n)
}

# Sigma = I_p.
cov_identity <- function(p) {
  diag(p)
}

# Compound-symmetric: Sigma_ii = 1, Sigma_ij = rho.
cov_compound_symmetric <- function(p, rho) {
  if (!is.numeric(rho) || length(rho) != 1) stop("rho must be a scalar.")
  m <- matrix(rho, nrow = p, ncol = p)
  diag(m) <- 1
  m
}

# AR(1): Sigma_ij = rho^|i-j|.
cov_ar1 <- function(p, rho) {
  if (!is.numeric(rho) || length(rho) != 1) stop("rho must be a scalar.")
  rho <- as.numeric(rho)
  stats::toeplitz(rho ^ (0:(p - 1)))
}

# Invert diagonal covariance, replacing zeros with tiny_constant.
invert_diagonal <- function(diag_values, tiny_constant = 1e-2) {
  d <- as.numeric(diag_values)
  d[d == 0] <- tiny_constant
  diag(1 / d)
}

# Diagonal LDA (uses only feature-wise variances).
dlda_fit <- function(x, y, tiny_constant = 1e-2) {
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("x and y have incompatible sizes.")

  p <- ncol(x)
  mu <- matrix(0, nrow = p, ncol = 2)

  if (p == 1) {
    mu[, 1] <- mean(x[y == 1, , drop = TRUE])
    mu[, 2] <- mean(x[y == -1, , drop = TRUE])
    s_diag <- diag(stats::cov(x))
  } else {
    mu[, 1] <- colMeans(x[y == 1, , drop = FALSE])
    mu[, 2] <- colMeans(x[y == -1, , drop = FALSE])
    s_diag <- diag(stats::cov(x))
  }

  s_inv <- invert_diagonal(s_diag, tiny_constant = tiny_constant)
  dif_vec <- matrix(mu[, 1] - mu[, 2], ncol = 1)
  sum_vec <- matrix(mu[, 1] + mu[, 2], ncol = 1)

  coef <- s_inv %*% dif_vec
  thresh <- 0.5 * t(coef) %*% sum_vec

  list(coef = coef, thresh = thresh, prec_matrix = s_inv, group_mean = mu)
}

# Predict {-1, 1} from a dlda_fit model.
predict_dlda <- function(newdata, model) {
  newdata <- as.matrix(newdata)
  score <- as.vector(newdata %*% model$coef)
  ifelse(score > as.numeric(model$thresh), 1, -1)
}

# Theoretical LDA error rate given population and sample statistics.
# mu1, mu2: population mean vectors (length p).
# sigma: population covariance matrix (p x p).
# xbar1, xbar2: sample mean vectors (length p).
# s_inv: inverse (or pseudoinverse) of the pooled sample covariance (p x p).
lda_error_rate <- function(mu1, mu2, sigma, xbar1, xbar2, s_inv) {
  mid <- 0.5 * (xbar1 + xbar2)
  dif <- xbar1 - xbar2
  sd_inv_dif <- as.vector(s_inv %*% dif)
  denom <- sqrt(as.numeric(t(dif) %*% s_inv %*% sigma %*% sd_inv_dif))
  if (denom < .Machine$double.eps) return(0.5)
  num1 <- -as.numeric(t(mu1 - mid) %*% sd_inv_dif)
  num2 <-  as.numeric(t(mu2 - mid) %*% sd_inv_dif)
  0.5 * stats::pnorm(num1 / denom) + 0.5 * stats::pnorm(num2 / denom)
}

# Wang (2018) error formula (equal class sizes), gamma = p/n in (0,1].
wang_error_equal_n <- function(delta_sq, gamma) {
  stats::pnorm(-(sqrt(1 - gamma) * (delta_sq / (2 * sqrt(delta_sq + 4 * gamma)))))
}

# Theorem 1 visualization helper. maha_scaled = Delta^2/n, gamma = p/n.
theorem1_viz <- function(maha_scaled, gamma) {
  stats::pnorm(-(2 * sqrt(maha_scaled) * sqrt(gamma * (1 - gamma))))
}
