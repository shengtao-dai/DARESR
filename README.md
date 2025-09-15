# DARESR

Distributionally Adaptive Regression Estimator under Series Expansion
## Installation

```r
# Option 1: pak (recommended)
install.packages("pak")
pak::pak("shengtao-dai/DARESR")

# Option 2: remotes
install.packages("remotes")
remotes::install_github("shengtao-dai/DARESR")

# Option 3: devtools
install.packages("devtools")
devtools::install_github("shengtao-dai/DARESR")

```

## Usage

```r
library(DARESR)

#Suppose Y (numeric), X (matrix/data.frame), and x0 (numeric vector)
res <- series_eq_estimate(Y, X, x0, p = 4, degree = 4, c_band = 1.0, K_grid = c(1,4,9), L = 5)
str(res)
```

