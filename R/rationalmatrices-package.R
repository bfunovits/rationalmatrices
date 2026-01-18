# ===================================================================
# Functional Area #10: Utilities & Helpers (CLAUDE.md)
#
# Purpose: Package-level documentation and initialization
#   - roxygen2 documentation directives
#   - Rcpp/RcppArmadillo integration
#   - Package imports and exports
#
# Related Files: All other files (package configuration)
# ===================================================================

#' @useDynLib rationalmatrices
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

# rationalmatrices.R

#' A Collection of Tools for Rational Matrices
#'
#' The package \code{rationalmatrices} collects classes, methods and tools 
#' for creating and manipulating rational matrices, i.e. matrices whose 
#' entries are rational functions. 
#' 
#' @section Classes:
#' 
#' There exist many different representations for rational matrices. This package in particular 
#' deals with "left matrix fraction" representations (implemented as \code{\link{lmfd}} class) and 
#' "statespace" representations (implemented as \code{\link{stsp}} class). As a special case of course 
#' also polynomial matrices (\code{\link{polm}} class) are implemented. The coefficients 
#' of the power series expansion of the rational function are stored as \code{\link{pseries}} objects and 
#' a collection of values of the rational function may be stored as  \code{\link{zvalues}} objects. 
#' 
#' The package in particular offers tools to convert one representation into another 
#' (equivalent) representation, see e.g. \code{\link{pseries2stsp}} and \code{\link{pseries2lmfd}}. 
#' 
#' @section Methods:
#' 
#' The package provides "standard" matrix - tools: 
#'
#' \itemize{
#' \item{General methods, like \code{\link{print}}, \code{\link{dim}}, \code{\link{plot}}, ...}
#' \item{Arithmetic operations, like addition and multiplication. See in particular 
#'      \code{\link[rationalmatrices]{Ops.ratm}} and \code{\link{\%r\%}}.}
#' \item{Extract parts of the matrix, transposition, (column or row) bind two or more matrices, ...} 
#' }
#' 
#' Some specific methods for rational matrices are 
#' 
#' \itemize{
#' \item{Compute poles and zeroes of rational matrices.}
#' \item{Check properties of the rational matrix, like stability and inverse stability.}
#' \item{Some support for the echelon canonical form, see e.g. \link{Kronecker-Indices}.}
#' \item{Normal Forms for polynomial matrices, like the Hermite normal form and the Smith form. 
#'       See \code{\link{hnf}}, \code{\link{snf}}, \code{\link{whf}} and \code{\link{col_reduce}}. 
#'       Check for "left (co)primeness" with \code{\link{is.coprime}}}. 
#' \item{Compute the derivative (with respect to \eqn{z}).}
#' \item{Tools related to statespace representations, e.g. computation of 
#'       controllability and observability Grammians and the computation of balanced 
#'       (truncated) realizations. See \code{\link{grammians}}, \code{\link{balance}}, ...}
#' }       
#'       
#' @section Author(s):
#' 
#' Wolfgang Scherrer and Bernd Funovits
#'
#' Maintainer: <bernd.funovits@gmail.com>
#'
#' References: See Also: Examples:
#' 
       
#' @keywords internal
"_PACKAGE"

#' @useDynLib rationalmatrices, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rdpack reprompt
NULL#' Pipe operator
#'
#' Re-export pipe operator \code{\%>\%} to turn function composition into a series of imperative statements.
#' For more extensive description, see function \code{`\%>\%`} in package \emph{magrittr}.
#'
#' @importFrom magrittr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs First argument of the function of the right-hand-side of the pipe operator.
#' @param rhs Function whose first argument is given by the left-hand-side argument \code{lhs} of the pipe operator.
#'
#' @examples
#' x = array(stats::rnorm(2*1*3, sd = 0.01), dim = c(2,1,3))
#' # Instead of
#' pseries(polm(x))
#' # you can write
#' x %>% polm() %>% pseries()
NULL
