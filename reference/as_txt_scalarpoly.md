# Coerce Scalar Polynomials to Character Strings

This utility coerces a scalar polynomial (given by the vector of real
coefficients) or a scalar Laurent polynomial to a character string. If
the vector should be interpreted as a Laurent polynomial, the minimal
degree should be given in the argument `laurent` described below. The
following "formats" are implemented.

- `syntax = "txt"` returns a simple text representation,

- `syntax = "TeX"` renders the coefficients to string in "TeX" syntax
  and

- `syntax = "expression"` gives a string which may be rendered to an `R`
  expression with [`parse`](https://rdrr.io/r/base/parse.html). This
  expression may be used to evaluate the polynomial and for annotating
  plots, see [`plotmath`](https://rdrr.io/r/grDevices/plotmath.html) and
  the examples below.

## Usage

``` r
as_txt_scalarpoly(
  coefs,
  syntax = c("txt", "TeX", "expression"),
  x = "z",
  laurent = FALSE
)
```

## Arguments

- coefs:

  Vector of doubles (complex elements are not allowed) representing a
  univariate polynomial (first element corresponds to power zero). When
  it represents a Laurent polynomial, the first element corresponds to
  the minimal degree.

- syntax:

  (character string) determines the format of the output string.

- x:

  (character string) polynomial variable.

- laurent:

  Boolean or integer. Default set to FALSE. If one deals with Laurent
  polynomials, an integer corresponding to the minimal degree of the
  Laurent polynomial should be supplied.

## Value

Character string used for printing univariate (Laurent) polynomials.

## Examples

``` r
coefs = c(1, 2.3, 0, -1, 0)

as_txt_scalarpoly(coefs, syntax = 'txt', x = 'x')
#> [1] "1 + 2.3x - x^3"
as_txt_scalarpoly(coefs, syntax = 'TeX', x = '\\alpha')
#> [1] "1 + 2.3\\alpha - \\alpha^{3}"
as_txt_scalarpoly(coefs, syntax = 'expression', x = 'z')
#> [1] "1 + 2.3*z - z^{3}"
as_txt_scalarpoly(coefs = sample((-10):10, 7, replace = TRUE), 
                  syntax = 'txt', x = 'x', 
                  laurent = -3)
#> [1] "8x^-3 - x^-2 + 2x^-1 + 6 + 3x + 2x^2"

if (FALSE) { # \dontrun{
# the case syntax = "expression" may be used e.g. as follows

# make_polyfun creates a "closure" which evaluates the polynomial at given points
# note that this simple version does not work for zero polynomials!
make_polyfun = function(coefs) {
  expr = parse(text = as_txt_scalarpoly(coefs, 'expression', 'x'))
  fun = function(x) {
    return(eval(expr))
  }
  return(fun)
}

a = make_polyfun(coefs)
a(1)   # return the value  of the polynomial at x = 1
a(1:5) # return the values of the polynomial at x = 1,2,3,4,5

# create a plot
x_grid = seq(from = -1, to = 1, length.out = 101)
plot(x_grid, a(x_grid), type = 'l', xlab = 'x', ylab = 'a(x)',
     main = parse(text = paste('a(x) == ',
          as_txt_scalarpoly(coefs, syntax = 'expression', x = 'x'))))
} # }
```
