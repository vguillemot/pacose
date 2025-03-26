#' The iPACOSE algorithm for graph estimation
#'
#' The iPACOSE algorithm allowing to iteratively estimate a graph from a dataset and a first estimation of the partial correlation matrix.
#'
#' @usage ipacose(x = x, pc = pc, method = "pacose.ridge", cutoff, gr = NULL, k = 2, cv.method = "CV", adaptive = FALSE)
#'
#' @param x a dataset (matrix) of dimensions n x p.
#' @param pc a first version of the p x p partial correlation matrix.
#' @param method should be "pacose.pls", "pacose.ridge" or "pacose.adalasso".
#' @param cutoff a threshold to apply to the new partial correlation matrices obtained with PACOSE in order to obtain a new version of the graph.
#' @param gr the true graph (an object of class [igraph::igraph]) when it is known. It must be set to NULL (the default value) when the true graph is not known.
#' @param k integer, the number of folds to be used when selecting the optimal value for the regularization parameter.
#' @param cv.method Only valid when the reference method is [ridge.net]: determines the way the ridge parameter is computed ("HKB" or "CV", see \code{pacose.ridge}).
#' @param adaptive boolean, default to FALSE.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Niter}{number of iterations.}
#'   \item{gr_it}{the final estimated graph.}
#'   \item{pcor_it}{the final estimated p x p partial correlation matrix.}
#'   \item{sens}{when the true graph is known, the sensitivity of the estimated graph.}
#'   \item{ppv}{when the true graph is known, the PPV of the estimated graph.}
#' }
#'
#' @references Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.
#'
#' @author Vincent Guillemot
#'
#' @examples
#' require(mvtnorm)
#' require(igraph)
#'
#' omega <- matrix(c(1     , -0.477, 0.304, 0.478, -0.591, 0,
#'                   -0.477, 2     , 0.206, 0    , 0.382 , 0,
#'                   0.304 , 0.206 , 1    , 0    , 0.181 , 0.242,
#'                   0.478 , 0     , 0    , 3    , 0.141 , 0,
#'                   -0.591, 0.382 , 0.181, 0.141, 1     , 0,
#'                   0     , 0     , 0.242, 0    , 0     , 2),
#'                 nrow = 6, ncol = 6)
#' gr <- graph.adjacency((omega != 0), mode = "undirected", diag = FALSE)
#' x <- rmvnorm(n = 100, sigma = solve(omega))
#'
#' # First estimation with ridge.net, the threshold is, for this toy example,
#' # arbitrarily set to 0.05
#' pcor0.0 <- ridge.net(x)$pcor
#' pcor0 <- (abs(pcor0.0) > 0.05) + 0
#'
#' # Use iPACOSE to estimate iteratively the partial correlation matrix
#' pcor1 <- ipacose(x = x, pc = pcor0, method = "pacose.ridge",
#'                  cutoff = 0.05, gr = NULL, cv.method = "HKB")
#' cat("Number of iterations: ", pcor1$Niter)
#'
#' invcov2pcor(omega)
#' pcor1$pcor_it
#'
#' @keywords algebra multivariate

ipacose <- function(x=x,pc=pc,method="pacose.ridge",cutoff,gr=NULL,k=2, cv.method="CV",adaptive=FALSE) {
  nonstop <- TRUE
  p <- ncol(x)
  iter <- 0
  gr1 <- gr2 <- graph.adjacency((abs(pc) >= cutoff) - diag(p),mode="undirected")
  while(nonstop) {
   iter <- iter+1
   res2 <- do.call(method,list(X=x,gg=gr2,k=k, cv.method=cv.method))
   if (method == "pacose.adalasso") {
      if (adaptive==TRUE){ pcor <- res2$pcor.lasso
      }else{ pcor <- res2$pcor.adalasso }
   }else{
      pcor <- res2$pcor
   }
   #if (is.null(cutoff) & iter==1) {
   #  cutoff <- fdrtool(pc[upper.tri(pc)], statistic = "correlation", plot = F,verbose=F)$param[1]
   #  print("Estimated cutoff = " ,cutoff,"\n")
   #}
   # print(round(c(cutoff,nbri(gr,gr2)/nbr(gr),nbri(gr,gr2)/nbr(gr2)),3))
   if (cutoff==0) {
      grnew <- gr2
   }else{
      grnew <- graph.adjacency( (abs(pcor) >= cutoff) - diag(p) ,mode="undirected")
   }
   nonstop <- nbri(grnew,gr2) != nbru(grnew,gr2)
   gr2 <- grnew
  }
  if (is.null(gr)) {
    sens <- ppv <- NULL
  } else {
    sens <- nbri(gr,gr2)/nbr(gr); ppv <- nbri(gr,gr2)/nbr(gr2)
  }
  return(list(Niter=iter,gr_it=gr2,pcor_it=pcor,sens=sens,ppv=ppv))
}

