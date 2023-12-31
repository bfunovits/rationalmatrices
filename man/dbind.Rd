% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tools.R
\name{dbind}
\alias{dbind}
\title{Bind Arrays}
\usage{
dbind(d = 1, ...)
}
\arguments{
\item{d}{Integer. Concatenate arrays along the dimension \code{d}.}

\item{...}{arrays.}
}
\value{
Array
}
\description{
\code{dbind(d, x, y, ...)} concatenates/binds an arbitrary number of arrays \code{x, y,...} along
the dimension \code{d}.
For matrices, \code{dbind(d = 1, x, y, ...)} is (essentially) equivalent to
\code{\link{rbind}} and \code{dbind(d = 2, x, y, ...)} corresponds to \code{\link{cbind}}.
If the number of dimensions of an argument, \code{x} say,  is less than \code{d}, then this argument
is treated as an array of dimension \code{c(dim(x),1,..,1)}.
}
\details{
The procedure makes some effort to keep the \code{dimnames} attribute of the arguments.
}
\examples{
x = test_array(dim = c(2,3,1))
y = test_array(dim = c(2,3,1), dimnames = TRUE)
z = test_array(dim = c(2,3,3), dimnames = TRUE)

# Bind along dimension 1 (row-binding for matrices)
dbind(1, x)
dbind(1, x, y)

# Bind along dimension 2 (col-binding for matrices)
dbind(2, x, y)

# Bind along dimension 3
dbind(3, x, y)
dbind(3, x, y, z)

# Example that throws an error
\dontrun{
dbind(1, x, y, z) # throws an error, since the array x,y,z are not compatible
}
}
