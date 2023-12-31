% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tools.R
\name{set_default_par}
\alias{set_default_par}
\title{Set Default Graphical Parameters}
\usage{
set_default_par(m = 1, n = 1)
}
\arguments{
\item{m, n}{integers. The procedure assumes that the device is
split up into an m-by-n array of (sub-) figures and scales
\code{cex.axis}, \code{cex.lab} and \code{cex.main}
correspondingy.}
}
\value{
invisible list with the original graphical parameters.
}
\description{
Set default graphical parameters.
}
\examples{
opar = set_default_par(2,3)
print(opar)
}
\keyword{internal}
