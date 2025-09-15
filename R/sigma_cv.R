construct_sigma_matrix <- function(Y, X_bs, Q_hat_mat, tau_grid, x_bs, h) {
  n <- nrow(X_bs)
  K <- length(tau_grid)
  M <- crossprod(X_bs) / n
  x_norm2 <- sum(x_bs^2)
  
  U_inv_list <- lapply(seq_len(K), function(j) {
    residuals_j <- abs(Y - Q_hat_mat[, j])
    mask <- residuals_j <= h
    X_selected <- X_bs[mask, , drop = FALSE]
    if (nrow(X_selected) == 0) {
      U_j <- diag(1e-1, ncol(X_bs))
    } else {
      U_j <- crossprod(X_selected) / (2 * n * h)
      if (is.nan(base::rcond(U_j)) || base::rcond(U_j) < 1e-12) {
        U_j <- U_j + diag(1e-6, ncol(U_j))
      }
    }
    solve(U_j)
  })
  
  tau_mat        <- outer(tau_grid, tau_grid, "*")
  kernel_matrix  <- outer(tau_grid, tau_grid, pmin) - tau_mat
  
  x_col <- matrix(as.numeric(x_bs), ncol = 1)
  V <- do.call(cbind, lapply(U_inv_list, function(U_inv) drop(U_inv %*% x_col)))
  
  quad_matrix <- t(V) %*% M %*% V
  Sigma <- kernel_matrix * quad_matrix / x_norm2
  Sigma
}


construct_sigma_matrix_l1 <- function(Y, X_bs, Q_hat_mat, tau_grid, x_bs, h) {
  n <- nrow(X_bs)
  K <- length(tau_grid)
  M <- crossprod(X_bs) / n
  x_norm2 <- sum(x_bs^2)

  U_inv_list <- lapply(seq_len(K), function(j) {
    residuals_j <- abs(Y - Q_hat_mat[, j])
    mask <- residuals_j <= h
    X_selected <- X_bs[mask, , drop = FALSE]
    if (nrow(X_selected) < ncol(X_bs)) {
      U_j <- diag(1e-6, ncol(X_bs))
    } else {
      U_j <- crossprod(X_selected) / (2 * n * h)
      if (is.nan(rcond(U_j)) || rcond(U_j) < 1e-12) {
        U_j <- U_j + diag(1e-6, ncol(U_j))
      }
    }
    solve(U_j)
  })

  tau_mat <- outer(tau_grid, tau_grid, "*")
  kernel_matrix <- outer(tau_grid, tau_grid, pmin) - tau_mat
  V <- sapply(U_inv_list, function(U_inv) U_inv %*% t(x_bs))
  quad_matrix <- t(V) %*% M %*% V
  Sigma <- kernel_matrix * quad_matrix / x_norm2
  Sigma
}

cv_select_K_series <- function(Y, X, K_grid, L = 5, h, p, degree) {
  if (is.vector(X)) X <- matrix(X, ncol = ncol(as.matrix(X)))
  n <- length(Y)
  set.seed(1L)
  folds <- sample(rep(seq_len(L), length.out = n))
  MSE_mat <- matrix(NA_real_, nrow = L, ncol = length(K_grid))

  for (l in seq_len(L)) {
    S_l <- which(folds == l)
    S_minus_l <- which(folds != l)

    X_train <- X[S_minus_l, , drop = FALSE]
    Y_train <- Y[S_minus_l]
    X_test  <- X[S_l, , drop = FALSE]
    Y_test  <- Y[S_l]

    knots_info_tr <- get_knots_info(X_train, p, degree)
    X_train_bs <- construct_bs_given_knots(X_train, p, degree, knots_info_tr)
    X_test_bs  <- construct_bs_given_knots(X_test,  p, degree, knots_info_tr)

    for (k_idx in seq_along(K_grid)) {
      K <- K_grid[k_idx]
      tau_grid <- seq_len(K) / (K + 1)

      coeff_hat_mat <- sapply(tau_grid, function(tau) series_quantile_coeff(Y_train, X_train_bs, tau))
      Q_hat_mat_test <- X_test_bs %*% coeff_hat_mat

      w_eq <- matrix(0, nrow = nrow(X_test_bs), ncol = K)
      for (i in seq_len(nrow(X_test_bs))) {
        x0_bs <- X_test_bs[i, , drop = FALSE]
        Sigma <- construct_sigma_matrix(Y_test, X_test_bs, Q_hat_mat_test, tau_grid, x0_bs, h)
        Sigma_reg <- Sigma + diag(1e-8, K)
        Sigma_inv <- solve(Sigma_reg)
        num_w <- rowSums(Sigma_inv)
        deno_w <- sum(Sigma_inv)
        w_eq[i, ] <- 0.5 * (num_w + rev(num_w)) / deno_w
      }

      m_hat <- rowSums(Q_hat_mat_test * w_eq)
      MSE_mat[l, k_idx] <- mean((Y_test - m_hat)^2)
    }
  }

  avg_MSE <- colMeans(MSE_mat)
  K_opt <- K_grid[which.min(avg_MSE)]
  list(K_opt = K_opt, avg_MSE = avg_MSE, MSE_mat = MSE_mat)
}

cv_select_K_series_l1 <- function(Y, X, K_grid, L = 5, h, p, degree) {
  if (is.vector(X)) X <- matrix(X, ncol = ncol(as.matrix(X)))
  n <- length(Y)
  set.seed(1L)
  folds <- sample(rep(seq_len(L), length.out = n))
  MSE_mat <- matrix(NA_real_, nrow = L, ncol = length(K_grid))

  for (l in seq_len(L)) {
    S_l <- which(folds == l)
    S_minus_l <- which(folds != l)

    X_train <- X[S_minus_l, , drop = FALSE]
    Y_train <- Y[S_minus_l]
    X_test  <- X[S_l, , drop = FALSE]
    Y_test  <- Y[S_l]

    knots_info_tr <- get_knots_info(X_train, p, degree)
    X_train_bs <- construct_bs_given_knots(X_train, p, degree, knots_info_tr)
    X_test_bs  <- construct_bs_given_knots(X_test,  p, degree, knots_info_tr)

    for (k_idx in seq_along(K_grid)) {
      K <- K_grid[k_idx]
      tau_grid <- seq_len(K) / (K + 1)

      coeff_hat_mat <- sapply(tau_grid, function(tau) series_quantile_coeff(Y_train, X_train_bs, tau))
      Q_hat_mat_test <- X_test_bs %*% coeff_hat_mat

      w_eq <- matrix(0, nrow = nrow(X_test_bs), ncol = K)
      for (i in seq_len(nrow(X_test_bs))) {
        x0_bs <- X_test_bs[i, , drop = FALSE]
        Sigma <- construct_sigma_matrix_l1(Y_test, X_test_bs, Q_hat_mat_test, tau_grid, x0_bs, h)
        Sigma_reg <- Sigma + diag(1e-8, K)
        Sigma_inv <- solve(Sigma_reg)
        num_w <- rowSums(Sigma_inv)
        deno_w <- sum(Sigma_inv)
        w_eq[i, ] <- 0.5 * (num_w + rev(num_w)) / deno_w
      }

      m_hat <- rowSums(Q_hat_mat_test * w_eq)
      MSE_mat[l, k_idx] <- mean((Y_test - m_hat)^2)
    }
  }

  avg_MSE <- colMeans(MSE_mat)
  K_opt <- K_grid[which.min(avg_MSE)]
  list(K_opt = K_opt, avg_MSE = avg_MSE, MSE_mat = MSE_mat)
}
