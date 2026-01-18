# Prune (Laurent) Matrix Polynomial

Performs three steps to simplify a matrix polynomial.

1.  All leading coefficients where the absolute values of the real and
    the imaginary parts are less than or equal to `tol` are set to zero.
    For Laurent polynomials, the same is happening from the other
    direction.

2.  The zero leading coefficient matrices are dropped. For Laurent
    polynomials, the same is happening for the coefficient matrices
    pertaining to low powers.

3.  If all the absolute values of the imaginary parts of the
    coefficients are less than or equal to `tol` then the coefficients
    are set to real values.

Empty polynomial matrices (i.e. matrices with zero rows or columns) are
set to polynomial of zero degree.

## Usage

``` r
prune(x, tol = sqrt(.Machine$double.eps), brutal = FALSE)
```

## Arguments

- x:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  or
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  object.

- tol:

  Double. Tolerance parameter. Default set to
  `sqrt(.Machine$double.eps)`.

- brutal:

  Boolean. Default set to FALSE. If TRUE, all small elements are set to
  zero (irrespective of whether they are leading coefficients or not).

## Value

A Laurent matrix polynomial, i.e. a
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
or
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object.

## Examples

``` r
x = polm(array(c(1,0,0,0,
                 0,1,0,0,
                 0,0,1,0,
                 0,0,0,1,
                 0,0,0,0), dim = c(2,2,5)) + 1e-4)
x
#> ( 2 x 2 ) matrix polynomial with degree <= 4 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]   [,2] z^3 [,1]   [,2] z^4 [,1]
#> [1,]   1.0001 1e-04   0.0001 1e-04    1e-04 1.0001    1e-04 0.0001    1e-04
#> [2,]   0.0001 1e-04   1.0001 1e-04    1e-04 0.0001    1e-04 1.0001    1e-04
#>       [,2]
#> [1,] 1e-04
#> [2,] 1e-04
prune(x, tol = 1e-3)
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]   [,2] z^3 [,1]   [,2]
#> [1,]   1.0001 1e-04   0.0000 1e-04        0 1.0001        0 0.0000
#> [2,]   0.0001 1e-04   1.0001 1e-04        0 0.0001        0 1.0001
prune(x, tol = 1e-3, brutal = TRUE)
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]   [,2] z^3 [,1]   [,2]
#> [1,]   1.0001     0   0.0000     0        0 1.0001        0 0.0000
#> [2,]   0.0000     0   1.0001     0        0 0.0000        0 1.0001

# Case of complex variables:
x = x + complex(imaginary = 1e-5)
x
#> ( 2 x 2 ) matrix polynomial with degree <= 4 
#>           z^0 [,1]         [,2]  z^1 [,1]     [,2] z^2 [,1]      [,2] z^3 [,1]
#> [1,] 1.0001+1e-05i 1e-04+1e-05i 0.0001+0i 1e-04+0i 1e-04+0i 1.0001+0i 1e-04+0i
#> [2,] 0.0001+1e-05i 1e-04+1e-05i 1.0001+0i 1e-04+0i 1e-04+0i 0.0001+0i 1e-04+0i
#>           [,2] z^4 [,1]     [,2]
#> [1,] 0.0001+0i 1e-04+0i 1e-04+0i
#> [2,] 1.0001+0i 1e-04+0i 1e-04+0i
prune(x, tol = 1e-3)
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]   [,2] z^3 [,1]   [,2]
#> [1,]   1.0001 1e-04   0.0000 1e-04        0 1.0001        0 0.0000
#> [2,]   0.0001 1e-04   1.0001 1e-04        0 0.0001        0 1.0001

# also works for constant matrix polynomials (i.e. matrices)
x = polm(array(0:3, dim = c(2,2,1))+1e-4)
x
#> ( 2 x 2 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]   [,2]
#> [1,]   0.0001 2.0001
#> [2,]   1.0001 3.0001
prune(x, tol = 1e-3)
#> ( 2 x 2 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]   [,2]
#> [1,]   0.0000 2.0001
#> [2,]   1.0001 3.0001

# empty polynomials are coerced to polynomials of degree zero 
x = polm(array(0, dim = c(0,2,5)))
x
#> ( 0 x 2 ) matrix polynomial with degree <= 4 
prune(x)
#> ( 0 x 2 ) matrix polynomial with degree <= -1 

# Laurent polynomials:
(lp = lpolm(c(0, 1:3, 0), min_deg = 2))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 6, and minimal degree >= 2
#>      z^2 [,1] z^3 [,1] z^4 [,1] z^5 [,1] z^6 [,1]
#> [1,]        0        1        2        3        0
prune(lp)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 5, and minimal degree >= 3
#>      z^3 [,1] z^4 [,1] z^5 [,1]
#> [1,]        1        2        3
```
