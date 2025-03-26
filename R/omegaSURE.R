#' SURE estimator of the covariance matrix under constraints
#'
#' The SURE of the covariance matrix under constraints is an estimate of the covariance matrix knowing that some elements of the precision matrix are equal to 0. These elements are given in the form of a graph.
#'
#' @usage omegaSURE(x, gr)
#'
#' @param x matrix, the dataset.
#' @param gr \code{\link[igraph:igraph-package]{igraph}} object representing the zeros in the inverse covariance matrix.
#'
#' @return \code{KhatSURE}: the estimate of the precision matrix knowing a pattern of 0s.
#'
#' @references Wiesel, A., Eldar, Y. C., and Hero, A. O. (2010). Covariance estimation in decomposable Gaussian graphical models. \emph{IEEE Transactions on Signal Processing}, 58(3):1482--1492.
#'
#' @author Vincent Guillemot
#'
#' @seealso \code{\link{omegaMVUE}}
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
#' gr <- simplify(graph.adjacency((omega != 0), mode = "undirected", diag = FALSE))
#' x <- rmvnorm(n = 100, sigma = solve(omega))
#' omega.hat <- omegaSURE(x, gr)
#'
#' omega
#' round(omega.hat, 3)
#'
#' @keywords algebra multivariate

omegaSURE <- function(x, gr) {
  if (!(is.chordal(gr,fillin=TRUE)$chordal)) {warning("Non decomposable graph");return(matrix(NaN,ncol(x),ncol(x)))}
  m <- get.adjacency(gr,sparse=F) ; colnames(m) <- rownames(m) <- 1:ncol(m)
  jt <- ug.to.jtree(m)
  l1 <- lapply(jt@cliques,function(u) as.numeric(u$vset))
  l2 <- lapply(jt@separators,function(u) as.numeric(u$separator))
  p <- ncol(x)
  normD2 <- sum((Reduce("+",lapply(l1,oneK,x)) - Reduce("+",lapply(l2,oneK,x)))**2)
  num <-  -2*( sum(sapply(l1,trk,x)) - sum(sapply(l2,trk,x)) )
  d <- num/normD2
  KhatSURE <- Reduce("+",lapply(l1,oneK_SURE,x,d)) - Reduce("+",lapply(l2,oneK_SURE,x,d))
  return(KhatSURE=KhatSURE)
}

