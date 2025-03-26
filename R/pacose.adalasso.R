#' The LASSO and adaptive LASSO versions of the PACOSE algorithm
#'
#' The function takes as an argument a dataset and a graph and returns an estimation of the partial correlation matrix.
#'
#' @usage pacose.adalasso(X, gg, k = 10, use.Gram = FALSE, both = TRUE, verbose = FALSE, cv.method = "CV")
#'
#' @param X a dataset (matrix) of dimensions n x p.
#' @param gg the graph (an object of class \code{\link[igraph:igraph-package]{igraph}}) to integrate.
#' @param k integer, the number of folds to be used when selecting the optimal value of the regularization parameter.
#' @param use.Gram see \code{\link[parcor:adalasso]{adalasso}}.
#' @param both boolean, whether to compute both the adaptive LASSO and the non adaptive LASSO versions of the partial correlation matrix or not, default to TRUE.
#' @param verbose boolean, whether to print out intermediate messages or not, default to FALSE.
#' @param cv.method equals "CV", for a cross validation determination of the regularization parameter. Alternative values are only used in function \code{\link{pacose.ridge}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{pcor.lasso}{an estimate of the p x p partial correlation matrix with the non adaptive version of the method LASSO.}
#'   \item{invcov.lasso}{an estimate of the p x p precision matrix with the non adaptive version of the method LASSO.}
#'   \item{pcor.adalasso}{an estimate of the p x p partial correlation matrix with the adaptive version of the method LASSO.}
#'   \item{invcov.adalasso}{an estimate of the p x p precision matrix with the adaptive version of the method LASSO.}
#' }
#'
#' @references Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.
#'
#' @author Vincent Guillemot
#'
#' @seealso \code{\link[parcor:adalasso.net]{adalasso.net}}, \code{\link{pacose.ridge}}, \code{\link{pacose.pls}}
#'
#' @examples
#' \dontrun{
#' require(mvtnorm)
#' require(igraph)
#'
#' omega <- matrix(c(1     , -0.477, 0.304, 0.478, -0.591, 0,
#'                    -0.477, 2     , 0.206, 0    , 0.382 , 0,
#'                    0.304 , 0.206 , 1    , 0    , 0.181 , 0.242,
#'                    0.478 , 0     , 0    , 3    , 0.141 , 0,
#'                    -0.591, 0.382 , 0.181, 0.141, 1     , 0,
#'                    0     , 0     , 0.242, 0    , 0     , 2),
#'                  nrow = 6, ncol = 6)
#' gr <- graph.adjacency((omega != 0), mode = "undirected", diag = FALSE)
#' x <- rmvnorm(n = 100, sigma = solve(omega))
#' omega.hat <- pacose.adalasso(x, gr)$invcov.lasso
#'
#' omega
#' round(omega.hat, 3)
#' }
#'
#' @keywords algebra multivariate


pacose.adalasso <- function (X, gg, k = 10, use.Gram = FALSE, both = TRUE, verbose = FALSE, cv.method="CV")  {
    p <- ncol(X)
    X <- scale(X)
    colnames(X) <- 1:p
    A <- get.adjacency(gg)
    B.lasso <- B.adalasso <- matrix(0, nrow = p, ncol = p)
    colnames(B.lasso) <- colnames(B.adalasso) <- 1:p
    pcor.adalasso <- NULL
    pv <- pvar.shrink(X,verbose=F)
    if (verbose) cat(paste("Performing local (adaptive) lasso regressions\n"))
    if (verbose) cat(paste("Vertex no "))
    for (i in 1:p) {
        if ((i/10) == floor(i/10)) {
            if (verbose) cat(paste(i, "..."))
        }
        noti <- (1:p)[-which(A[i,] == 0)]
        #noti <- (1:p)[-i]
        yi <- X[, i]
        Xi <- X[, noti]
        if(length(noti)==0) { #pv[i] <- var(yi)
        } else {
          if (!is.null(dim(Xi))) {
            dummy <- adalasso(Xi, yi, k = k, use.Gram = use.Gram,both = both)
            coefi.lasso <- dummy$coefficients.lasso
            B.lasso[i, noti] <- coefi.lasso
            if (both == TRUE) {
               coefi.adalasso <- dummy$coefficients.adalasso
               B.adalasso[i, noti] <- coefi.adalasso
            }
          } else {
              lmo <- lm(yi~Xi)
              B.lasso[i, noti ] <- lmo$coefficients[2]
              if (both == TRUE) {
                 B.adalasso[i, noti] <- lmo$coefficients[2]
              }
              #pv[i] <- var(lmo$residuals)
          }

        }
    }
    pcor.lasso <- Beta2parcor(B.lasso, verbose = verbose)
    invcov.lasso <- pcor2invcov(pcor.lasso,pv)
    if (both == TRUE) {
        pcor.adalasso <- Beta2parcor(B.adalasso, verbose = verbose)
        invcov.adalasso <- pcor2invcov(pcor.adalasso,pv)
    }
    if (verbose) cat(paste("\n"))
    return(list(pcor.lasso = pcor.lasso, pcor.adalasso = pcor.adalasso,
                invcov.lasso = invcov.lasso, invcov.adalasso = invcov.adalasso))
}

