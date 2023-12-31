% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_methods.R
\name{as.lpolm}
\alias{as.lpolm}
\alias{as.lpolm.polm}
\title{Coerce to Laurent polynom object}
\usage{
as.lpolm(obj, ...)

\method{as.lpolm}{polm}(obj, ...)
}
\arguments{
\item{obj}{\code{\link{polm}} object}

\item{...}{other arguments}
}
\value{
\code{\link{lpolm}} object
}
\description{
The attribute \code{min_deg} is set to zero for the given function input.
}
\examples{
p = test_polm(degree = 2)
as.lpolm(p)
}
