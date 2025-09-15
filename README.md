# DARESR

Distributionally Adaptive Regression Estimator under Series Expansion
## Installation

```r
# Install devtools if needed
# install.packages("devtools")

devtools::install_local("DARESR", upgrade = "never")
# or build first:
# devtools::build("DARESR"); devtools::install("DARESR_0.1.0.tar.gz")
```

## Usage

```r
library(DARESR)

#Suppose Y (numeric), X (matrix/data.frame), and x0 (numeric vector)
res <- series_eq_estimate(Y, X, x0, p = 4, degree = 4, c_band = 1.0, K_grid = c(1,4,9), L = 5)
str(res)
```

