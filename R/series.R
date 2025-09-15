series_cond <- function(Y, X_bs) {
  fit <- stats::lm(Y ~ X_bs - 1)
  as.numeric(X_bs %*% stats::coef(fit))
}

series_cond_point <- function(Y, X_bs, x_bs) {
  x_bs <- matrix(x_bs, nrow = 1)
  fit <- stats::lm(Y ~ X_bs - 1)
  as.numeric(x_bs %*% stats::coef(fit))
}

series_quantile <- function(Y, X_bs, tau) {
  fit <- quantreg::rq(Y ~ X_bs - 1, tau = tau)
  as.numeric(X_bs %*% stats::coef(fit))
}

series_quantile_point <- function(Y, X_bs, x_bs, tau) {
  x_bs <- as.matrix(x_bs)
  fit <- quantreg::rq(Y ~ X_bs - 1, tau = tau)
  as.numeric(x_bs %*% stats::coef(fit))
}

series_quantile_coeff <- function(Y, X_bs, tau) {
  fit <- quantreg::rq(Y ~ X_bs - 1, tau = tau)
  stats::coef(fit)
}
