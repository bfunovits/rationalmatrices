% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_methods.R
\name{as.rmfd}
\alias{as.rmfd}
\alias{as.rmfd.pseries}
\title{Coerce to Right Matrix Fraction Description}
\usage{
as.rmfd(obj, method, ...)

\method{as.rmfd}{pseries}(obj, method, ...)
}
\arguments{
\item{obj}{object}

\item{method}{character string}

\item{...}{optional additional arguments}
}
\value{
object of class \code{\link{rmfd}}
}
\description{
The function \code{\link{as.rmfd.pseries}} calls \code{\link{pseries2rmfd}} with default parameters.
Of course the \code{\link{pseries}} object must contain sufficiently many lags.
}
