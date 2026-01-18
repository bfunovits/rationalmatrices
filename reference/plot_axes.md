# Plot Grid Lines, Axes and Axes Labels

This tools plots the grid lines,axes and axes titles/labels. The x-axis
may be a date/time axis.

## Usage

``` r
plot_axes(
  axes = c(TRUE, TRUE, FALSE, FALSE),
  tick_labels = axes,
  x_date = 0,
  titles = rep(NA_character_, 4),
  style = style_parameters("bw"),
  parse_titles = FALSE
)
```

## Arguments

- axes:

  4-dimensional boolean vector, indicating whether an axis should be put
  at the respective side of the plot.

- tick_labels:

  4-dimensional boolean vector, indicating whether the respective axis
  has tick-labels.

- x_date:

  scalar number, used to convert the (numeric) x-axis limits `xlim` to a
  date/time object. If `x_date` is of class `POSIXct` then the limits
  are converted as `as.POSIXct(xlim, origin = x_date)`. If `x_date` is
  of class `Date` then the limits are converted as
  `as.Date(xlim, origin = x_date)`. Otherwise no coercion is performed.
  Finally `graphics::Axis(x = xlim, ...)` is called to do the actual
  plotting of the x-axis.

- titles:

  a 4-dimensional character vector. If the respective entry is `NA` then
  no label/title is put at this side of the plot.

- style:

  list of "style" parameters as returned by
  [`style_parameters`](https://bfunovits.github.io/rationalmatrices/reference/style_parameters.md).

- parse_titles:

  if TRUE the titles are coerced to
  [`expression`](https://rdrr.io/r/base/expression.html) before
  plotting.

## Value

(invisible) `NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
date = seq(as.Date("2019/1/1"), by = "month", length.out = 24)
titles = c(NA, 'y','main', NA)
plot(date, 1:24, axes = FALSE, xlab = NA, ylab = NA)
plot_axes(axes = c(TRUE, FALSE, FALSE, TRUE), titles = titles,
          style = style_parameters('gray'), x_date = date[1] - as.numeric(date[1]))
points(date, 1:24)
} # }
```
