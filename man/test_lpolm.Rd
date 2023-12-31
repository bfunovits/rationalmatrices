% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{test_lpolm}
\alias{test_lpolm}
\title{Create Test Laurent Polynomial Matrix}
\usage{
test_lpolm(
  dim = c(1, 1),
  degree_max = 1,
  degree_min = -1,
  random = FALSE,
  col_start_matrix = NULL,
  value_at_0 = NULL,
  col_end_matrix = NULL
)
}
\arguments{
\item{dim}{two dimensional vector of non negative integers, determines the dimension of 
the polynomial matrix to be created. If the prescribed number of rows or number of columns 
is zero then an "empty" polynomial matrix is generated. In this case all 
parameters below are ignored.}

\item{degree_max, degree_min}{Integer, vector, or matrix specifying minimal and maximal (column) degrees}

\item{random}{If TRUE the coefficents are generated by drawing from a normal distribution.
If FALSE then the coefficient of the \eqn{k}-th power \eqn{z^k} of 
the (i,j)-th entry is set equal to "\eqn{ijk}". In this case the 
parameters below are ignored!}

\item{col_start_matrix}{Force a matrix for lowest column degrees}

\item{value_at_0}{desired value of the polynomial at \eqn{z=0} (or \code{NULL}).}

\item{col_end_matrix}{desired column end matrix (or \code{NULL}).}
}
\value{
\code{\link{lpolm}} object
}
\description{
See \code{\link{test_polm}}.
This simple tool creates Laurent polynomial matrices for testing purposes.
The \code{column_start_matrix}, i.e. the matrix consisting of the columns pertaining to the smallest column degrees,
the \code{value_at_0}, and that \code{column_end_matrix} matrix are set in this order.
Thus, \code{column_end_matrix} may potentially overwrite \code{value_at_0}.
}
\details{
In contrast to \code{\link{test_polm}}, the absolute value of the zeros cannot be forced to be larger than a certain lower bound \code{bzeroes}.
}
\examples{
(lp = test_lpolm(dim = c(1,1), degree_max = 1, degree_min = -2))
(lp = test_lpolm(dim = c(3,3), degree_max = c(0,1,2), degree_min = -2))
(lp = test_lpolm(dim = c(3,3), degree_max = 1, degree_min = c(0,-1,-2)))
}
