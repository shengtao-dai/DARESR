#' @export
series_eq_estimate <- function(Y, X, x0, p, degree, c_band = 1.0, K_grid, L = 5) {
  # 统一成矩阵，列数是维度 d
  if (is.vector(X)) X <- matrix(X, ncol = length(x0))
  
  n <- nrow(X)
  h <- c_band * n^(-1/5)
  
  # 1) 先用训练集 X 的分位点生成 knots/boundary
  knots_info <- get_knots_info(X, p, degree)
  
  # 2) 基于同一套 knots 构造 X 与 x0 的样条基（保证维度一致）
  X_bs  <- construct_bs_given_knots(X,                p, degree, knots_info)
  x0_bs <- construct_bs_given_knots(matrix(x0, 1L),   p, degree, knots_info)
  
  # 3) CV 选择 K
  cv_result <- cv_select_K_series(Y, X, K_grid, L, h, p, degree)
  K_opt <- cv_result$K_opt
  if (K_opt < 1L) stop("cv_select_K_series must return K_opt >= 1.")
  
  # 4) 计算 K 个分位点的拟合值（整样本）
  tau_grid  <- seq_len(K_opt) / (K_opt + 1)
  Q_hat_mat <- sapply(tau_grid, function(t) series_quantile(Y, X_bs, t))  # n × K
  
  # 5) 用 x0 的基构造 Sigma 并求权重
  Sigma     <- construct_sigma_matrix(Y, X_bs, Q_hat_mat, tau_grid, x0_bs, h)
  Sigma_reg <- Sigma + diag(1e-8, nrow = K_opt, ncol = K_opt)
  Sigma_inv <- solve(Sigma_reg)
  
  num_w  <- rowSums(Sigma_inv)
  deno_w <- sum(Sigma_inv)
  w_eq   <- 0.5 * (num_w + rev(num_w)) / deno_w
  
  # 6) 计算 x0 处每个分位的点估计，再加权
  Q_hat_mat_x <- sapply(tau_grid, function(t) series_quantile_point(Y, X_bs, x0_bs, t))
  m_hat_series_eq <- sum(Q_hat_mat_x * w_eq)
  
  # 7) 简单方差近似（与原代码保持一致）
  M2 <- sum(x0_bs^2)
  cov_hat            <- sqrt(n * deno_w / M2)
  cov_hat_series_eq  <- 1 / cov_hat
  m_hat_series_eq_up <- m_hat_series_eq + 1.96 / cov_hat_series_eq
  m_hat_series_eq_lo <- m_hat_series_eq - 1.96 / cov_hat_series_eq
  
  list(
    m_hat_series_eq        = m_hat_series_eq,
    cov_hat_series_eq      = cov_hat_series_eq,
    m_hat_series_eq_lower  = m_hat_series_eq_lo,
    m_hat_series_eq_upper  = m_hat_series_eq_up
  )
}
