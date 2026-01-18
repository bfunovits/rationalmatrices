# Start a new Plot

This tools starts a new plot and writes the optional "labels" `xlab`,
`ylab` and `main` into the (outer) margins of the plot. In addition an
(optional) legend is put into the right (outer) margin. The outer
margins of the plot are determined based on whether or not text or
legend is put into the respective margin.

## Usage

``` r
start_plot(xlab = NULL, ylab = NULL, main = NULL, legend_args = NULL)
```

## Arguments

- xlab, ylab, main:

  (optional) character or
  [`expression`](https://rdrr.io/r/base/expression.html). If `NULL` then
  the respective outer margin is set to zero.

- legend_args:

  (optional) list with
  [`legend`](https://rdrr.io/r/graphics/legend.html) arguments. Note
  that the legend is always put at the "right" side of the plot. If
  `NULL` then the right outer margin is set to zero.

## Value

(invisible) 4-dimensional vector with the outer margins `par("mar")`.

## See also

[`legend`](https://rdrr.io/r/graphics/legend.html) and
[`par`](https://rdrr.io/r/graphics/par.html).

## Examples

``` r
if (FALSE) { # \dontrun{
set_default_par()
start_plot(xlab = 'x', ylab = 'y', main = 'main',
           legend_args = list(legend = c('eins','zwei'), fill = c('red','blue')))
           graphics::box(which = 'inner', col = 'red')
} # }
```
