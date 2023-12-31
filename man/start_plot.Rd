% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tools.R
\name{start_plot}
\alias{start_plot}
\title{Start a new Plot}
\usage{
start_plot(xlab = NULL, ylab = NULL, main = NULL, legend_args = NULL)
}
\arguments{
\item{xlab, ylab, main}{(optional) character or \code{\link{expression}}.
If \code{NULL} then the respective outer margin is set to zero.}

\item{legend_args}{(optional) list with \code{\link[graphics]{legend}}
arguments. Note that the legend is always put at the "right" side of the
plot. If \code{NULL} then the right outer margin is set to zero.}
}
\value{
(invisible) 4-dimensional vector with the outer margins
        \code{par("mar")}.
}
\description{
This tools starts a new plot and writes the optional "labels"
\code{xlab}, \code{ylab} and \code{main} into the (outer) margins
of the plot. In addition an (optional) legend is put into the
right (outer) margin. The outer margins of the plot are determined
based on whether or not text or legend is put into the respective margin.
}
\examples{
\dontrun{
set_default_par()
start_plot(xlab = 'x', ylab = 'y', main = 'main',
           legend_args = list(legend = c('eins','zwei'), fill = c('red','blue')))
           graphics::box(which = 'inner', col = 'red')
}
}
\seealso{
\code{\link[graphics]{legend}} and \code{\link[graphics]{par}}.
}
\keyword{internal}
