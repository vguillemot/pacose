
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pacose

`pacose` estimates partial correlation and precision matrices when a
graph specifies conditional independence constraints. It provides Ridge,
PLS, LASSO, adaptive LASSO, iPACOSE, MVUE, SURE and INVEST estimators.

## Installation

Install the development version from GitHub:

``` r
remotes::install_github("vguillemot/pacose")
```

## Example

The graph below encodes the allowed non-zero entries of a precision
matrix. The Ridge version estimates the constrained precision matrix
from simulated observations.

``` r
library(igraph)
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union
library(pacose)

set.seed(1)
omega <- matrix(c(1, .4, 0, 0,
                  .4, 1, .3, 0,
                  0, .3, 1, .2,
                  0, 0, .2, 1), nrow = 4, byrow = TRUE)
graph <- graph_from_adjacency_matrix(omega != 0, mode = "undirected",
                                     diag = FALSE)
X <- mvtnorm::rmvnorm(100, sigma = solve(omega))

fit <- pacose.ridge(X, graph, cv.method = "HKB")
round(fit$invcov, 2)
#>      [,1] [,2] [,3] [,4]
#> [1,] 1.03 0.45 0.00 0.00
#> [2,] 0.45 1.00 0.36 0.00
#> [3,] 0.00 0.36 1.03 0.15
#> [4,] 0.00 0.00 0.15 1.05
```

## Reference

Guillemot V., Bender A., Boulesteix A.-L. (2013). *Iterative
Reconstruction of High-Dimensional Gaussian Graphical Models Based on a
New Method to Estimate Partial Correlations under Constraints*. **PLOS
ONE**, 8(4), e60536. <https://doi.org/10.1371/journal.pone.0060536>
