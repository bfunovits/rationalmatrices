# Plot 3D Arrays

This internal helper function plots data arranged in 3-dimensional
arrays. The base application is as follows. Let x,y be two 3-dimensional
arrays with dimension c(m,n,l). Then `plot_3D` splits the device region
into a a matrix like array of (m times n) subfigures. In the (i,j)-th
subfigure `y[i,j,]` is plotted against `x[i,j,]`.

## Usage

``` r
plot_3D(
  x,
  y,
  xlim = c("subfig", "column", "global"),
  ylim = c("subfig", "row", "global"),
  log = "",
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  subfigure_main = "(i_,j_)-th entry",
  parse_subfigure_main = FALSE,
  style = c("gray", "bw", "bw2", "colored"),
  col = NA,
  type = "l",
  lty = "solid",
  lwd = 1,
  pch = 16,
  cex.points = 1,
  bg.points = "black",
  legend = NULL,
  legend_args = NA
)
```

## Arguments

- x:

  list of 3D-arrays (define the 'x'-values)

- y:

  list of 3D-arrays (define the 'y'-values)

- xlim, ylim:

  determine the axis limits of the subfigures. E.g. `xlim = 'column'`
  means that all subfigures in a column use the same x-axis limits. The
  parameter `xlim` may also contain a 2-dimensional vector `c(x1,x2)`.
  In this case all sub-figures use the given limits for the x-axis.
  Furthermore the limits for the y-axis are computed based on the
  corresponding data "subset".

- log:

  a character string which contains "x" if the x axis is to be
  logarithmic, "y" if the y axis is to be logarithmic and "xy" or "yx"
  if both axes are to be logarithmic.

- main:

  (character or [`expression`](https://rdrr.io/r/base/expression.html))
  main title of the plot

- xlab:

  (character string or
  [`expression`](https://rdrr.io/r/base/expression.html)) label for the
  x-axis

- ylab:

  (character or [`expression`](https://rdrr.io/r/base/expression.html))
  label for the y-axis

- subfigure_main:

  scalar or `(m x n)` matrix of type "character" with the titles for the
  subfigures. If subfigure_main is a scalar character string then the
  procedures creates a matrix of respective titles by replacing the
  "place holders" '`i_`' and '`j_`' with the respective row and column
  number.

- parse_subfigure_main:

  boolean. If `TRUE` then the titles for the subfigures are parsed to
  [`expression`](https://rdrr.io/r/base/expression.html) before
  plotting. See also
  [`plotmath`](https://rdrr.io/r/grDevices/plotmath.html) on the usage
  of expressions for plot annotations.

- style:

  (character string) determines the appearance of the plot (background
  color of the plot regions, color and line style of the grid lines,
  axis color, ...) See also
  [`style_parameters`](https://bfunovits.github.io/rationalmatrices/reference/style_parameters.md).

- col:

  vector of line colors

- type:

  vector of plot types. The following values are possible: "p" for
  points, "l" for lines, "b" for both points and lines, "c" for empty
  points joined by lines, "o" for overplotted points and lines, "s" and
  "S" for stair steps and "h" for histogram-like vertical lines. `'n'`
  suppresses plotting.

- lty:

  vector of line types. Line types can either be specified as integers
  (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash,
  5=longdash, 6=twodash) or as one of the character strings "blank",
  "solid", "dashed", "dotted", "dotdash", "longdash", or "twodash",
  where "blank" uses ‘invisible lines’ (i.e., does not draw them).

- lwd:

  vector of line widths.

- pch:

  vector of plotting character or symbols. See
  [`points`](https://rdrr.io/r/graphics/points.html) for possible
  values.

- cex.points:

  vector of scales for the plotting symbols.

- bg.points:

  vector of fill color for the open plot symbols.

- legend:

  (character or [`expression`](https://rdrr.io/r/base/expression.html)
  vector). If `NULL` then no legend is produced.

- legend_args:

  (optional) list with parameters for the legend. A legend title can be
  included with `legend_args = list(title = my_legend_title)`. Note that
  the slots `x`, `y` are ignored and the legend is always put at the
  right hand side of the plot. See also
  [`legend`](https://rdrr.io/r/graphics/legend.html).

## Value

A [`function`](https://rdrr.io/r/base/function.html) ("closure"),
`subfig` say, which may be to add additional graphical elements to the
subfigures. The call `opar = subfig(i,j)` sets up the coordinate system
margins and figure coordinates such that one may add lines, text,
points, ... to the `(i,j)`-th subfigure.

## Details

This function is mainly used by the plot methods, see
[`plot methods`](https://bfunovits.github.io/rationalmatrices/reference/plot.md).

For more information about the parameters
`col, type, lty, ...,bg.points` see
[`plot.default`](https://rdrr.io/r/graphics/plot.default.html).

This helper function only makes some very basic checks on the given
parameters.
