#' Transforms a partial correlation matrix into an inverse covariance matrix
#'
#' Transforms an input p x p positive definite partial correlation matrix into an inverse covariance matrix given the vector of the partial variances.
#'
#' @usage pcor2invcov(A, pv)
#'
#' @param A a given partial correlation matrix.
#' @param pv a vector of partial variances.
#'
#' @return The inverse covariance matrix deduced from the given partial correlations and variances.
#'
#' @references Whittaker, J. (1990). *Graphical models in applied multivariate statistics*. Wiley.
#'
#' @keywords algebra multivariate
#' @export
pcor2invcov <- function(A, pv) {
   # if (any(is.nan(A))) return(NA)
   A <- -A ; diag(A) <- -diag(A)
   vec <- 1/pv
   sqrt(diag(vec))%*%A%*%sqrt(diag(vec))
}

