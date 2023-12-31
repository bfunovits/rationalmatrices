% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_tools.R
\name{btext}
\alias{btext}
\title{Write "Boxed" Text into the Margins of a Plot.}
\usage{
btext(
  text,
  side = 3,
  cex = 1,
  col = "black",
  size = cex,
  bg = "lightgray",
  border = "lightgray",
  parse_text = FALSE
)
}
\arguments{
\item{text}{character or \code{\link{expression}}.}

\item{side}{on which side of the plot (1=bottom, 2=left, 3=top, 4=right).}

\item{cex}{character expansion factor.}

\item{col}{color to use (for the text).}

\item{size}{determines the size of the box (in "line" units)}

\item{bg, border}{background and border color of the box.}

\item{parse_text}{if yes the procedure tries to coerce the
text into an \code{\link{expression}}.}
}
\value{
(invisible) vector with the "box" corner and center coordinates.
}
\description{
This tool writes text into a "colored" box in one of the margins
of the plot region.
}
\examples{
\dontrun{
plot(1:10, axes = FALSE)
graphics::box()
btext('hallo', side = 1)
btext('Sigma[11]', side = 2, bg = 'lightblue')
btext('Sigma[11]', side = 3, bg = NA, border = 'black', size = 1.2, parse_text = TRUE)
btext(expression(Sigma[11]), side = 4, bg = 'orange', border = 'black', cex = 1.5)
}
}
\seealso{
\code{\link[graphics]{mtext}} and \code{\link[grDevices]{plotmath}}
for details on mathematical annotation.
}
\keyword{internal}
