# Coerce Scalar Polynomial Filters to Character Strings

This utility coerces a scalar polynomial filter (given by the vector of
coefficients) to a character string. The following "formats" are
implemented. `syntax = "txt"` returns a simple text representation,
`syntax = "TeX"` renders the coefficients to string in "TeX" syntax and
`syntax = "expression"` gives a string which may be rendered to an `R`
expression with [`parse`](https://rdrr.io/r/base/parse.html). This
expression may be used to evaluate the filter and for annotating plots,
see [`plotmath`](https://rdrr.io/r/grDevices/plotmath.html) and the
examples below.

## Usage

``` r
as_txt_scalarfilter(
  coefs,
  syntax = c("txt", "TeX", "expression"),
  x = "z",
  t = "t"
)
```

## Arguments

- coefs:

  (numeric) vector of coefficients.

- syntax:

  (character string) determines the format of the output string.

- x:

  (character string) names the "input" series.

- t:

  (character string) names the "time-index".

## Value

character string.

## See also

[`as_txt_scalarpoly`](https://bfunovits.github.io/rationalmatrices/reference/as_txt_scalarpoly.md)
and
[`as_tex_matrixfilter`](https://bfunovits.github.io/rationalmatrices/reference/as_tex_matrixfilter.md).

## Examples

``` r
coefs = c(1, 2.3, 0, -1, 0)

as_txt_scalarfilter(coefs, syntax = 'txt', x = 'x', t = 't')
#> [1] "x[t] + 2.3x[t-1] - x[t-3]"
as_txt_scalarfilter(coefs, syntax = 'TeX', x = 'x', t = 's')
#> [1] "x_{s} + 2.3x_{s-1} - x_{s-3}"
as_txt_scalarfilter(coefs, syntax = 'expression', x = 'x', t = 'k')
#> [1] "x[k] + 2.3*x[k-1] - x[k-3]"

if (FALSE) { # \dontrun{
# the case syntax = "expression" may be used e.g. as follows

# make_filterfun creates a "closure" which computes the filter-output
# note that this simple version does not work for zero filters!
make_filterfun = function(coefs) {
  p = length(coefs) - 1
  expr = parse(text = as_txt_scalarfilter(coefs, 'expression', 'x', 't'))
  fun = function(x, t) {
    # x, t must be vectors
    y = rep(NA_real_, length(t))
    t0 = t
    y = rep(NA_real_, length(t))

    i = ((t0 > p) & (t0 <= length(x)))
    t = t0[i]
    if (any(i)) y[i] = eval(expr)

    return(y)
  }
  return(fun)
}

coefs = rep(1, 4) / 4  # represents a moving average of length 4.
a = make_filterfun(coefs)
u = rnorm(100)       # input series
a(u, 1)    # return the value of the output series at t = 1
           # this value is not defined due to missing initial values
a(u, 1:10) # return the values of the output series at t = 1,..,10

# create a plot
plot(1:length(u), u, type = 'n', xlab = 'time', ylab = '')
grid()
lines(1:length(u), u, col = 'black', lwd = 1)
lines(1:length(u), a(u, 1:length(u)), col = 'red', lwd = 2)
legend('topright', bty = 'n',
       fill = c('black', 'red'),
       legend = c(expression(u[t]),
                  parse(text = paste('x[t] == ',
                     as_txt_scalarfilter(coefs, 'expression', 'u','t')))) )
} # }
```
