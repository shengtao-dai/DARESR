# Distributionally adaptive regression estimator (series expansion)

Estimates the conditional mean \\m(x_0) = E\[Y \mid X = x_0\]\\ by a
weighted combination of series quantile regression estimates at the
levels \\\tau_k = k / (K + 1)\\, \\k = 1, \dots, K\\. Each covariate is
expanded in a B-spline basis, conditional quantiles are fitted with
[`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html), and the
weights come from the inverse of the estimated asymptotic covariance
matrix of the K quantile estimates. K is chosen from `K_grid` by
`L`-fold cross-validation.

## Usage

``` r
series_eq_estimate(Y, X, x0, p, degree, c_band = 1, K_grid, L = 5)
```

## Arguments

- Y:

  Numeric response vector of length n.

- X:

  Numeric matrix or data frame of covariates (n x d).

- x0:

  Numeric vector of length d, the target point.

- p:

  Integer, number of equally spaced knot intervals per covariate (p - 1
  interior knots).

- degree:

  Integer, degree of the B-spline basis.

- c_band:

  Numeric scalar, the bandwidth constant.

- K_grid:

  Integer vector of candidate numbers of quantiles K.

- L:

  Integer, number of cross-validation folds.

## Value

A list with

- `m_hat_series_eq`: the estimate of \\m(x_0)\\;

- `cov_hat_series_eq`: its estimated standard error;

- `m_hat_series_eq_lower`, `m_hat_series_eq_upper`: the 95% confidence
  interval, `m_hat_series_eq` -/+ 1.96 \* `cov_hat_series_eq`.

## Details

The bandwidth used for the covariance estimate is \\h = c\_{band} \\
n^{-1/5}\\. Cross-validation folds are drawn after `set.seed(1)`, so
results are reproducible, but the global random number generator is
reset.

## References

Dai, S. (2025). Distributionally Adaptive Regression Estimators. Working
paper. <https://ssrn.com/abstract=5683602>

## Examples

``` r
set.seed(123)
n <- 500
X <- cbind(rnorm(n), rnorm(n))
Y <- as.numeric(0.5 + X %*% c(1, -1) + rt(n, df = 3))  # heavy-tailed errors
res <- series_eq_estimate(Y, X, x0 = c(0, 0), p = 4, degree = 4,
                          K_grid = c(1, 3), L = 2)
#> Warning: some 'x' values beyond boundary knots may cause ill-conditioned bases
#> Warning: some 'x' values beyond boundary knots may cause ill-conditioned bases
#> Warning: some 'x' values beyond boundary knots may cause ill-conditioned bases
res$m_hat_series_eq  # true value is 0.5
#> [1] 0.3585145
```
