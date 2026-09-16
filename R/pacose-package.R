#' @keywords internal
#' @useDynLib pacose
#' @importFrom corpcor pvar.shrink pseudoinverse
#' @importFrom glmnet glmnet
#' @importFrom igraph E graph.adjacency graph.intersection graph.union get.adjacency is.chordal
#' @importFrom lcd ug.to.jtree
#' @importFrom MASS lm.ridge
#' @importFrom ppls penalized.pls.cv
#' @importFrom stats coef cov cov2cor lm predict var
"_PACKAGE"

#' PACOSE and iPACOSE implementation plus several other covariance selection methods
#'
#' PACOSE is a method able to estimate a partial correlation matrix knowing a pattern of zeros. iPACOSE is an iterative application of PACOSE aiming at estimating a graph from a dataset.
#'
#' @name pacose-package
#'
#' @details
#' This R package contains the core functions presented by Guillemot et al. [1], along with some reference covariance selection methods that were not yet implemented in R. PACOSE is designed to estimate partial correlation matrices under constraints, and iPACOSE estimates a Gaussian graphical model [2] from a dataset.
#'
#' The following covariance selection methods are implemented in this package along with PACOSE and iPACOSE: \code{\link{omegaMVUE}}, \code{\link{omegaSURE}} [3], \code{\link{INVEST_wrapper}} [4].
#'
#' To illustrate our method, three different versions of PACOSE are implemented: \code{\link{pacose.ridge}}, \code{\link{pacose.pls}}, \code{\link{pacose.adalasso}}. They are based respectively on the Ridge, PLS and LASSO regressions.
#'
#' iPACOSE is implemented in the function \code{\link{ipacose}}.
#'
#' The functions \code{mylars}, \code{adalasso} and \code{ridge.cv} are internal local-regression helpers. The public PACOSE wrappers expose the corresponding estimation methods.
#'
#' @author
#' Vincent Guillemot \email{contact@vincentguillemot.fr}, Andreas Bender.
#'
#' @references
#' [1] Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.
#'
#' [2] Whittaker, J. (1990). Graphical models in applied multivariate statistics. Wiley.
#'
#' [3] Wiesel, A., Eldar, Y. C., and Hero, A. O. (2010). Covariance estimation in decomposable gaussian graphical models. IEEE Transactions on Signal Processing, 58(3):1482--1492.
#'
#' [4] Wermuth, N. and Scheidt, E. (1977). Fitting a covariance selection model to a matrix, algorithm 105. Journal of the Royal Statistical Society C, 26:88--92.
#'
#' [5] Kraemer, N., Schaefer, J., and Boulesteix, A.-L. (2009). Regularized estimation of large scale gene association networks using gaussian graphical models. BMC Bioinformatics, 10:384.
#'
#' @keywords package
#'
#' @seealso The public estimation functions are \code{\link{pacose.ridge}}, \code{\link{pacose.pls}} and \code{\link{pacose.adalasso}}.
#'
#' @examples
#' # For further examples on simulated and real datasets, see the associated
#' # article [Guillemot et al., 2012] and the companion website.
NULL
