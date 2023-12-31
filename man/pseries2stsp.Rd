% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_methods.R
\name{pseries2stsp}
\alias{pseries2stsp}
\alias{Ho-Kalman}
\title{Ho-Kalman Realization Algorithm}
\usage{
pseries2stsp(
  obj,
  method = c("balanced", "echelon"),
  Hsize = NULL,
  s = NULL,
  nu = NULL,
  tol = sqrt(.Machine$double.eps),
  Wrow = NULL,
  Wcol = NULL
)
}
\arguments{
\item{obj}{\code{\link{pseries}} object or 3-D array with dimension \eqn{(m,n,l+1)}.}

\item{method}{Character string, which determines the method and the "parametrization" type
of the state space model. See below for more details.}

\item{Hsize}{integer vector \code{c(f,p)}, number of block rows and block columns
of the Hankel matrix which is used to construct the statespace realization.
If NULL a default choice is made.}

\item{s}{desired state dimension. Only used for \code{method = "balanced"}.
Note however, if \eqn{s} is larger than the rank of the Hankel matrix,
then the procedure will break down.
If \code{s} is missing, then the state dimension is determined from
the singular values of the Hankel matrix.
To be precise the state dimension is chosen as the number of singular values which are
greater than or equal to \code{tol} times the maximum singular value.}

\item{nu}{Kronecker indices. Only used for  \code{method = "echelon"}.
If missing, then \code{nu} is computed with a QR decomposition of the
transpose of the Hankel matrix of the impulse response coefficients.}

\item{tol}{tolerance parameter used for the QR decomposition or the SVD decomposition
of the Hankel matrix \eqn{H} of the impulse response coefficients.}

\item{Wrow, Wcol}{weighting matrices (default is no weighting, i.e. identity matrices).
These weighting matrices are only used for \code{method="balanced"}, where the
SVD of the weighted Hankel matrix \code{Wrow \%*\% H \%*\% t(Wcol)} is computed.}
}
\value{
List with slots
\item{Xs}{\code{\link{stsp}} object, the rational matrix in statespace form}
\item{Hsv}{Singular values of the Hankel matrix for \code{method='balanced'} and \code{NULL} else.}
\item{nu}{Kronecker indices for \code{method='echelon'} and \code{NULL} else.}
}
\description{
This helper function implements the Ho-Kalman algorithm.
}
\details{
The procedure(s) may be used for model reduction (with some care).

There are a number of restrictions on the number of lags \eqn{l}
of the impulse response, the number of block rows (\eqn{f}),
block columns (\eqn{p}) of the Hankel matrix and the Kronecker
indices \eqn{\nu_i}{\nu[i]}. We require that:
\eqn{p>0}, \eqn{f>1},
\eqn{l \geq f+p-1}{l\ge f+p-1} and \eqn{\nu_i <f}{\nu[i] < f}.
If these restrictions are not satisfied an error is thrown.
}
\examples{
# generate random rational matrix X(z) in statespace form
# make sure that the A matrix is stable
m = 3
n = 2
s = 7
A = matrix(rnorm(s*s), nrow = s, ncol = s)
A = A / (1.1 * max(abs(eigen(A, only.values = TRUE)$values)))
Xs = stsp(A, B = matrix(rnorm(s*n), nrow = s, ncol = n),
          C = matrix(rnorm(s*m), nrow = m, ncol = s),
          D = diag(1, nrow = m, ncol = n))
Xi = pseries(Xs, lag.max = 20)

out = pseries2stsp(Xi, method = 'balanced')
print(out)
# check impulse response
all.equal(pseries(out$Xs, lag.max = 20), Xi)

Xs1 = as.stsp(Xi)
all.equal(Xs1, out$Xs)

out = pseries2stsp(Xi, method = 'echelon')
print(out)
# check impulse response
all.equal(pseries(out$Xs, lag.max = 20), Xi)

Xs1 = as.stsp(Xi, method = 'echelon')
all.equal(Xs1, out$Xs)
}
