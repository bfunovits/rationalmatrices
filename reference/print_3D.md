# Helper for Printing Matrix (Laurent) Polynomials

Called by the `rationalmatrices` methods for the
[print](https://rdrr.io/r/base/print.html) S3 generic.

## Usage

``` r
print_3D(
  a,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  laurent = FALSE
)
```

## Arguments

- a:

  Array with 3 dimensions, representing a polynomial or Laurent
  polynomial

- digits:

  Default NULL. Otherwise an integer specifying the number of digits of
  the coefficients that should be printed, see
  [`round`](https://rdrr.io/r/base/Round.html).

- format:

  Default "i\|jz", i.e. the matrix polynomial is printed that the "z
  part is the slowest moving index. The vertical bar designates
  separation across array dimensions (one bar = two array dimensions,
  two bars = array of dimension 3). In addition to multiple ordering
  options regarding row-index `i`, column-index `j`, and polynomial
  power `z`, there is also the option `character` (partial matching is
  enabled such that `format = "c"` gives the desired result) which
  prints each univariate (Laurent) polynomial into a matrix with
  appropriate numbers of rows and columns. Note that
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  objects have no format option. The option `'character'` is only
  implemented for polynomials, Laurent polynomials, LMFDs and RMFDs with
  real coefficients.

- laurent:

  Boolean or integer. Default set to FALSE. If one deals with Laurent
  polynomials, an integer corresponding to the minimal degree of the
  Laurent polynomial should be supplied.

## Value

Printed (Laurent) polynomial matrix
