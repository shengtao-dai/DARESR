#' Distributionally adaptive regression estimator (series expansion)
#'
#' Estimates the conditional mean \eqn{m(x_0) = E[Y \mid X = x_0]} by a
#' weighted combination of series quantile regression estimates at the levels
#' \eqn{\tau_k = k / (K + 1)}, \eqn{k = 1, \dots, K}. Each covariate is
#' expanded in a B-spline basis, conditional quantiles are fitted with
#' [quantreg::rq()], and the weights come from the inverse of the estimated
#' asymptotic covariance matrix of the K quantile estimates. K is chosen from
#' `K_grid` by `L`-fold cross-validation.
#'
#' The bandwidth used for the covariance estimate is
#' \eqn{h = c_{band} \, n^{-1/5}}. Cross-validation folds are drawn after
#' `set.seed(1)`, so results are reproducible, but the global random number
#' generator is reset.
#'
#' @param Y Numeric response vector of length n.
#' @param X Numeric matrix or data frame of covariates (n x d).
#' @param x0 Numeric vector of length d, the target point.
#' @param p Integer, number of equally spaced knot intervals per covariate
#'   (p - 1 interior knots).
#' @param degree Integer, degree of the B-spline basis.
#' @param c_band Numeric scalar, the bandwidth constant.
#' @param K_grid Integer vector of candidate numbers of quantiles K.
#' @param L Integer, number of cross-validation folds.
#' @return A list with
#'   * `m_hat_series_eq`: the estimate of \eqn{m(x_0)};
#'   * `cov_hat_series_eq`: its estimated standard error;
#'   * `m_hat_series_eq_lower`, `m_hat_series_eq_upper`: the 95% confidence
#'     interval, `m_hat_series_eq` -/+ 1.96 * `cov_hat_series_eq`.
#' @references Dai, S. (2025). Distributionally Adaptive Regression
#'   Estimators. Working paper. \url{https://ssrn.com/abstract=5683602}
#' @examples
#' set.seed(123)
#' n <- 500
#' X <- cbind(rnorm(n), rnorm(n))
#' Y <- as.numeric(0.5 + X %*% c(1, -1) + rt(n, df = 3))  # heavy-tailed errors
#' res <- series_eq_estimate(Y, X, x0 = c(0, 0), p = 4, degree = 4,
#'                           K_grid = c(1, 3), L = 2)
#' res$m_hat_series_eq  # true value is 0.5
#' @export
series_eq_estimate <- function(Y, X, x0, p, degree, c_band = 1.0, K_grid, L = 5) {
  if (is.vector(X)) X <- matrix(X, ncol = length(x0))
  
  n <- nrow(X)
  h <- c_band * n^(-1/5)
  
  knots_info <- get_knots_info(X, p, degree)
  
  X_bs  <- construct_bs_given_knots(X,                p, degree, knots_info)
  x0_bs <- construct_bs_given_knots(matrix(x0, 1L),   p, degree, knots_info)
  
  cv_result <- cv_select_K_series(Y, X, K_grid, L, h, p, degree)
  K_opt <- cv_result$K_opt
  if (K_opt < 1L) stop("cv_select_K_series must return K_opt >= 1.")
  
  tau_grid  <- seq_len(K_opt) / (K_opt + 1)
  Q_hat_mat <- sapply(tau_grid, function(t) series_quantile(Y, X_bs, t))  # n × K
  
  Sigma     <- construct_sigma_matrix(Y, X_bs, Q_hat_mat, tau_grid, x0_bs, h)
  Sigma_reg <- Sigma + diag(1e-8, nrow = K_opt, ncol = K_opt)
  Sigma_inv <- solve(Sigma_reg)
  
  num_w  <- rowSums(Sigma_inv)
  deno_w <- sum(Sigma_inv)
  w_eq   <- 0.5 * (num_w + rev(num_w)) / deno_w
  
  Q_hat_mat_x <- sapply(tau_grid, function(t) series_quantile_point(Y, X_bs, x0_bs, t))
  m_hat_series_eq <- sum(Q_hat_mat_x * w_eq)
  
  M2 <- sum(x0_bs^2)
  cov_hat            <- sqrt(n * deno_w / M2)
  cov_hat_series_eq  <- 1 / cov_hat
  m_hat_series_eq_up <- m_hat_series_eq + 1.96 * cov_hat_series_eq
  m_hat_series_eq_lo <- m_hat_series_eq - 1.96 * cov_hat_series_eq
  
  list(
    m_hat_series_eq        = m_hat_series_eq,
    cov_hat_series_eq      = cov_hat_series_eq,
    m_hat_series_eq_lower  = m_hat_series_eq_lo,
    m_hat_series_eq_upper  = m_hat_series_eq_up
  )
}
