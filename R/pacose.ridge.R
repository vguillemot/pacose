#' The Ridge version of the PACOSE algorithm
#'
#' The function takes as an argument a dataset and a graph and returns an estimation of the partial correlation matrix.
#'
#' @usage pacose.ridge(X, gg, lambda = NULL, scale = FALSE, k = 10, verbose = FALSE, cv.method = "CV")
#'
#' @param X a dataset (matrix) of dimensions n x p.
#' @param gg the graph (an object of class \code{\link[igraph:igraph-package]{igraph}}) to integrate.
#' @param lambda a set (vector) of values among which to choose when estimating the partial correlations.
#' @param scale boolean, whether to scale the data or not, default to FALSE.
#' @param k integer, the number of folds to be used when selecting the optimal value among \code{lambda}.
#' @param verbose boolean, whether to print out intermediate messages or not, default to FALSE.
#' @param cv.method determines the way the ridge parameter is computed: should be equal to either "HKB" for an analytical determination, or "CV", for a cross validation alternative. When equal to "HKB", neither \code{lambda} nor \code{k} are used.
#'
#' @return A list containing:
#' \describe{
#'   \item{pcor}{an estimate of the p x p partial correlation matrix.}
#'   \item{invcov}{an estimate of the p x p precision matrix.}
#' }
#'
#' @references Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.
#'
#' @author Vincent Guillemot
#'
#' @seealso \code{\link[parcor:ridge.net]{ridge.net}}, \code{\link{pacose.pls}}, \code{\link{pacose.adalasso}}
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
#' gr <- graph.adjacency(omega != 0, mode = "undirected", diag = FALSE)
#' x <- rmvnorm(n = 100, sigma = solve(omega))
#' omega.hat <- pacose.ridge(x, gr, cv.method = "HKB")$invcov
#'
#' omega
#' round(omega.hat, 3)
#'
#' @keywords algebra multivariate

pacose.ridge <- function(X, gg, lambda = NULL, scale = FALSE, k = 10, verbose = FALSE, cv.method="CV") {
    if (is.null(lambda) == TRUE) {
        ss <- seq(-10, -1, length = 1000)
        ss <- 10^ss
        n <- nrow(X)
        nn <- n - floor(n/k)
        lambda <- ss * nn * ncol(X)
    }
    n <- nrow(X)
    p <- ncol(X)
    pv <- pvar.shrink(X,verbose=F)
    X <- scale(X, scale = scale)
    A <- get.adjacency(gg)
    B <- matrix(0, nrow = p, ncol = p)
    lambda.opt <- rep(0, p)
    #pv <- rep(0, p)
    if (verbose) cat(paste("Performing local ridge regressions\n"))
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
              if (cv.method == "CV" ) lambda.opt.i <- ridge.cv(Xi, yi, lambda = lambda, scale = scale, plot.it = FALSE, k = k)$lambda.opt
              if (cv.method == "HKB") lambda.opt.i <- lm.ridge(yi~Xi, lambda = 0, scale = scale)$kHKB

              rr <- lm.ridge(yi ~ Xi, scale = scale, lambda = lambda.opt.i)

              B[i, noti ] <- coef(rr)[-1]
              lambda.opt[i] <- lambda.opt.i
              #pv[i] <- var(yi-(Xi%*%coef(rr)[-1]+coef(rr)[1]))
          } else {
              lmo <- lm(yi~Xi)
              B[i, noti ] <- lmo$coefficients[2]
              lambda.opt[i] <- 0
              #pv[i] <- var(lmo$residuals)
          }
        }
    }
    if (verbose) cat("\n")
    # print(lambda.opt)
    pcor <- Beta2parcor(B, verbose = verbose)
    invcov <- pcor2invcov(pcor,pv)
    return(list(pcor = pcor,invcov=invcov))
}

