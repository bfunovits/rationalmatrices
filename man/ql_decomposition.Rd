% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{QL and LQ decomposition}
\alias{QL and LQ decomposition}
\alias{ql_decomposition}
\alias{lq_decomposition}
\title{QL and LQ Decomposition}
\usage{
ql_decomposition(x, ...)

lq_decomposition(x, ...)
}
\arguments{
\item{x}{Matrix.}

\item{...}{Other arguments for \code{\link{qr}}}
}
\value{
List with two elements:
\item{q}{Semi-orthogonal matrix}
\item{l}{Lower triangular matrix with non-negative diagonal elements.}
}
\description{
Returns the QL and LQ decomposition of a matrix with non-negative "diagonal" elements in the QL and LQ decompositions.
Only works if no column pivoting occurs in the QR decomposition \code{\link{qr}}.
}
\note{
Shouldn't export this function when publishing package!
}
\section{Implementation of the QL decomposition using the QR decomposition}{

The base function \link{qr} is used in the following way.
(We need to assume that there is no pivoting since otherwise the function throws an error).
First, we reorder the columns of the input matrix \code{x} from last to first QR-decompose this matrix.
Next, we reorder the columns of \code{Q} and the rows of \code{R} in the same way. 
The original matrix \code{x} is now in \code{QL} decomposition.
}

\examples{
set.seed(1803)

# Tall matrix
x = matrix(stats::rnorm(5*3), 5, 3)
out = ql_decomposition(x)
all.equal(x, out$q \%*\% out$l)
all.equal(diag(ncol(out$q)), t(out$q) \%*\% out$q)
out$l
out = lq_decomposition(x)
all.equal(x, out$l \%*\% out$q)
all.equal(diag(nrow(out$q)), out$q \%*\% t(out$q))
out$l

# Wide matrix
x = matrix(stats::rnorm(5*3), 3, 5)
out = ql_decomposition(x)
all.equal(x, out$q \%*\% out$l)
all.equal(diag(ncol(out$q)), t(out$q) \%*\% out$q)
out$l
out = lq_decomposition(x)
all.equal(x, out$l \%*\% out$q)
all.equal(diag(nrow(out$q)), out$q \%*\% t(out$q))
out$l

# Square matrix
x = matrix(stats::rnorm(4*4), 4, 4)
out = ql_decomposition(x)
all.equal(x, out$q \%*\% out$l)
all.equal(diag(ncol(out$q)), t(out$q) \%*\% out$q)
out$l
out = lq_decomposition(x)
all.equal(x, out$l \%*\% out$q)
all.equal(diag(nrow(out$q)), out$q \%*\% t(out$q))
out$l

# reset seed
set.seed(NULL)
}
