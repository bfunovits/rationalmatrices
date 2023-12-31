% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polm_methods.R
\name{col_end_matrix}
\alias{col_end_matrix}
\title{Column End Matrix of a Polynomial Matrix}
\usage{
col_end_matrix(x)
}
\arguments{
\item{x}{A polynomial matrix, i.e. an object of class \code{\link{polm}}.}
}
\value{
The column end matrix.
}
\description{
The \emph{column end matrix} of an \eqn{(m,n)}-dimensional polynomial matrix 
\eqn{a(z)=a_0 + a_1 z + \cdots + a_p z^p}{a(z)=a[0] + a[1] z + \dots + a[p] z^p} is defined as follows. 
Suppose that the maximum degree of the elements in the \eqn{i}-th column is \eqn{p_i}{p[i]}. Then 
the column end matrix is the \eqn{(m,n)} matrix with \eqn{i}-th column equal to the 
\eqn{i}-th column of the coefficient matrix \eqn{a_{p_i}}{a[p[i]]}. If a column of 
\eqn{a(z)} is zero, then
the elements of the corresponding column of the column end matrix are set to \code{NA}'s.
}
\examples{
x = polm(array(c(0,1,1,0,
                 0,0,1,0,
                 0,0,0,1,
                 0,0,0,0), dim = c(2,2,4)))
x
degree(x)
degree(x, 'columns')
col_end_matrix(x)
}
