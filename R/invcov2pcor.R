#' Transforms an inverse covariance matrix into a partial correlation matrix
#'
#' Transforms the positive definite \eqn{p \times p} input precision matrix into the
#' corresponding \eqn{p \times p} partial correlation matrix.
#'
#' @usage invcov2pcor(mat)
#' @param mat the given inverse covariance matrix.
#' @return The corresponding partial correlation matrix.
#' @references Whittaker, J. (1990). *Graphical models in applied multivariate statistics*. Wiley.
#' @keywords algebra multivariate

invcov2pcor <- function(mat) cov2cor(oppdiag(mat))

