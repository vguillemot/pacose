#' Wrapper for the C implementation of Wermuth's covariance selection algorithm
#'
#' An R wrapper for a C implementation of Wermuth's covariance selection
#' algorithm. Takes as an input a dataset and the indices where there should
#' be zeros in the precision matrix and returns the corresponding precision
#' matrix (if possible) according to Wermuth's algorithm.
#'
#' @usage INVEST_wrapper(x, ind, delta = 1e-09, itermax = 2000)
#'
#' @param x a dataset (matrix) of dimensions n x p.
#' @param ind matrix of the indices corresponding to the zeros in the inverse
#'   covariance matrix.
#' @param delta positive real, value under which a coefficient in the inverse
#' covariance matrix is considered to be equal to 0.
#' @param itermax integer, a maximum number of iterations to avoid infinite
#'   loops.
#'
#' @return The desired inverse covariance matrix with zeros where specified.
#' @references Wermuth, N. and Scheidt, E. (1977).
#' Fitting a covariance selection model to a matrix, algorithm 105.
#' *Journal of the Royal Statistical Society C, 26*, 88--92.
#' @author Andreas Bender, Vincent Guillemot
#' @keywords algebra multivariate

INVEST_wrapper <- function(x, ind, delta = 10e-10, itermax = 2000) {
 tempmat <- matrix(10,ncol(x),ncol(x))
 maxiter <- itermax
 indices <- ind
 icovx <- pseudoinverse(cov(x))
 p <- ncol(x)
 invcovx_vecform <- .C("wermuthC",
   p = as.integer(p), ninteract=as.integer(nrow(indices)),
   delta = as.double(delta), error = as.double(1),
   iter = as.integer(0), maxiter = as.integer(maxiter),
   result = as.double(tempmat),PACKAGE="pacose")$result
 invcovx <-  matrix(invcovx_vecform, p, p)
 error   <- max(abs(invcovx[indices + 1]))
 return(list(invcovx = invcovx,error = error))
}

