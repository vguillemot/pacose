#' The PLS version of the PACOSE algorithm
#'
#' The function takes as an argument a dataset and a graph and returns an estimation of the partial correlation matrix.
#'
#' @usage pacose.pls(X, gg, scale = TRUE, k = 10, Ncomp = NULL, verbose = FALSE, cv.method = "CV")
#'
#' @param X a dataset (matrix) of dimensions n x p.
#' @param gg the graph (an object of class \code{\link[igraph:igraph-package]{igraph}}) to integrate.
#' @param scale boolean, whether to scale the data or not, default to TRUE.
#' @param k integer, the number of folds to be used when selecting the optimal value among \code{Ncomp}.
#' @param Ncomp integer, the maximal number of PLS components to be used.
#' @param verbose boolean, whether to print out intermediate messages or not, default to FALSE.
#' @param cv.method equals "CV", for a cross validation determination of the regularization parameter. Alternative values are only used in function \code{\link{pacose.ridge}}.
#'
#' @return A list containing:
#' \describe{
#'   \item{pcor}{an estimate of the p x p partial correlation matrix.}
#'   \item{invcov}{an estimate of the p x p precision matrix.}
#'   \item{m}{vector containing all the optimal values for the regularization parameters.}
#' }
#'
#' @references Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.
#'
#' @author Vincent Guillemot
#'
#' @seealso \code{\link[parcor:pls.net]{pls.net}}, \code{\link{pacose.ridge}}, \code{\link{pacose.adalasso}}
#'
#' @examples
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
#' omega.hat <- pacose.pls(x, gr, verbose = TRUE)$invcov
#'
#' omega
#' round(omega.hat, 3)
#'
#' @keywords algebra multivariate

pacose.pls <- function(X, gg, scale = TRUE, k = 10, Ncomp = NULL, verbose = FALSE, cv.method="CV") {
    n <- nrow(X)
    p <- ncol(X)
    pv <- pvar.shrink(X,verbose=F)
    k <- max(1,floor(k))
    if (k > n) {
        cat(paste("k exceeds the number of observations. Leave-one-out is applied.\n"))
        k <- n
    }
    A <- get.adjacency(gg)
    B <- matrix(0, nrow = p, ncol = p)
    m <- vector(length = p)
    if (verbose) cat(paste("Performing local pls regressions\n"))
    kernel <- FALSE
    if (n < (p - 1)) {
        kernel <- TRUE
    }
    if (verbose) cat(paste("Vertex no "))
    for (i in 1:p) {
        if ((i/10) == floor(i/10)) {
            if (verbose) cat(paste(i, "..."))
        }
        noti <- (1:p)[-which(A[i,] == 0)]
        yi <- X[, i]
        if(length(noti)==0) { #pv[i] <- var(yi)
        } else {
          Xi <- matrix(X[, noti],ncol=length(noti))
          Xi <- X[, noti]
          if (!is.null(dim(Xi))) {
              if (is.null(Ncomp)) {
                     ncomp <- min(n - 1, ncol(Xi))
              } else {
                     ncomp <- Ncomp
              }
              fit <- penalized.pls.cv(Xi, yi, scale = scale, k = k, ncomp = ncomp)
              B[i, noti] <- fit$coefficients
              m[i] <- fit$ncomp.opt
          } else {
              lmo <- lm(yi~Xi)
              B[i, noti ] <- lmo$coefficients[2]
              m[i] <- 0
              #pv[i] <- var(lmo$residuals)
          }
        }
    }
    if (verbose) cat("\n")
    pcor <- Beta2parcor(B, verbose = verbose)
    invcov <- pcor2invcov(pcor,pv)
    return(list(pcor = pcor,invcov=invcov, m = m))
}

