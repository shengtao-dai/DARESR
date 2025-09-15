test_that("basic run on small dat", {
  set.seed(123)
  n <- 500
  X <- cbind(rnorm(n), rnorm(n))
  beta <- c(1, -1)
  Y <- as.numeric(0.5 + X %*% beta + rnorm(n))
  x0 <- c(0, 0)
  K_grid <- c(1, 3)
  
  expect_error(series_eq_estimate(Y, X, x0, p = 4, degree = 4, c_band = 1.0, K_grid = K_grid, L=2), NA)
})
