construct_bs_basis <- function(X, p, degree) {
  delta <- 10^(-320)
  if (is.vector(X)) X <- matrix(X, nrow = 1)   # ★ 关键修复：1×d，而不是 d×1
  n <- nrow(X); d <- ncol(X)
  f <- p - 1 + degree
  X_bs <- matrix(0, nrow = n, ncol = d * f)
  for (j in 1:d) {
    xj <- X[, j]
    bounds <- stats::quantile(xj, probs = c(delta, 1 - delta), na.rm = TRUE)
    knots  <- seq(bounds[1], bounds[2], length.out = p + 1)[-c(1, p + 1)]
    bs_j <- splines::bs(xj, degree = degree, knots = knots,
                        Boundary.knots = bounds, intercept = FALSE)
    X_bs[, ((j - 1) * f + 1):(j * f)] <- bs_j
  }
  cbind(1, X_bs)
}


get_knots_info <- function(X, p, degree) {
  if (is.vector(X)) X <- matrix(X, nrow = 1)   # 容错
  d <- ncol(X)
  delta <- 10^(-320)
  knots_list <- vector("list", d)
  boundary_list <- vector("list", d)
  for (j in 1:d) {
    xj <- X[, j]
    bounds <- stats::quantile(xj, probs = c(delta, 1 - delta), na.rm = TRUE)
    knots  <- seq(bounds[1], bounds[2], length.out = p + 1)[-c(1, p + 1)]
    knots_list[[j]]    <- knots
    boundary_list[[j]] <- bounds
  }
  list(knots = knots_list, boundary = boundary_list)
}


construct_bs_given_knots <- function(X, p, degree, knots_info) {
  if (is.vector(X)) X <- matrix(X, nrow = 1)   # ★ 关键修复
  n <- nrow(X); d <- ncol(X)
  f <- p - 1 + degree
  X_bs <- matrix(0, nrow = n, ncol = d * f)
  for (j in 1:d) {
    xj <- X[, j]
    bs_j <- splines::bs(xj, degree = degree,
                        knots = knots_info$knots[[j]],
                        Boundary.knots = knots_info$boundary[[j]],
                        intercept = FALSE)
    X_bs[, ((j - 1) * f + 1):(j * f)] <- bs_j
  }
  cbind(1, X_bs)
}

