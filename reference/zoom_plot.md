# Zoom and Scroll

The utility `zoom_plot` creates and runs a "shiny app" which displays
two plots: Drawing (and dragging) a "brush" in the "control" plot at the
bottom zooms into zooms into the selected x-range of the plot at the
top.

## Usage

``` r
zoom_plot(p, p0 = p, title = NULL, ...)
```

## Arguments

- p, p0:

  two [`function`](https://rdrr.io/r/base/function.html) ("closure")
  objects which produce an R plot when calling `p()`, `p(xlim)` and
  `p0()` respectively. Both plots should use the same range of
  "x-values".

- title:

  (character) optimal title for the "shiny app window"

- ...:

  not used.

## Value

This function normally does not return; interrupt R to stop the
application (usually by pressing Ctrl+C or Esc).

## Details

This utility is based on the shiny package.

## Examples

``` r
if (FALSE) { # \dontrun{
make_plot_fun = function(x, y, type, col) {
   fun = function(xlim = NULL) {
   plot(x, y, xlim= xlim, type = type, col = col)
   }
   return(fun)
}

p = make_plot_fun(1:10, rnorm(10), type = 'p', col = 'red')
p0 = make_plot_fun(1:10, sin(1:10), type = 'l', col = 'black')
zoom_plot(p, p0, title = 'test "zoom_plot"')
} # }
```
