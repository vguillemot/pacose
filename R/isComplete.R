#' Test if a given (sub)graph is complete
#'
#' Determines if a given graph is complete or not.
#'
#' @usage isComplete(G, set)
#'
#' @param G an \code{\link[igraph:igraph-package]{igraph}} object.
#' @param set vector of the subgraph's nodes.
#'
#' @return Boolean: TRUE if the (sub)graph is complete, FALSE if not.
#'
#' @keywords graphs


isComplete <- function (G, set){
 if (length(set)<=1) {
   Res <- TRUE
 } else {
   B       <- as(G, "matrix")
   B       <- B[set,set]
   diag(B) <- 1
   Res     <- all(B != 0)
 }
 return(Res)
}
