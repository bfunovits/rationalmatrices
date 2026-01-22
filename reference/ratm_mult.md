# Matrix Multiplication of Rational Matrices

Matrix Multiplication of Rational Matrices

## Usage

``` r
e1 %r% e2
```

## Arguments

- e1, e2:

  two rational matrix objects (i.e.
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  objects), or objects which may be coerced to rational matrix objects.

## Value

rational matrix object. The class depends on the classes of the
arguments `e1,e1`.

## Examples

``` r
(lp1 = test_lpolm(degree_max = 1, degree_min = -1, random = TRUE))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1]  z^0 [,1] z^1 [,1]
#> [1,]  0.552598 0.2997174 1.005277
(lp2 = test_lpolm(degree_max = 1, degree_min = -1, random = TRUE))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>       z^-1 [,1]  z^0 [,1]   z^1 [,1]
#> [1,] -0.7531662 0.5506115 0.05285033
lp1 %r% lp2
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 2, and minimal degree >= -2
#>       z^-2 [,1]  z^-1 [,1]  z^0 [,1]  z^1 [,1]   z^2 [,1]
#> [1,] -0.4161981 0.07852988 -0.562908 0.5693574 0.05312924
```
