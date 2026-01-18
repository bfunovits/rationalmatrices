# Inflate a scalar to an (m x n) matrix with identical entries

Used in elementwise group operations
[`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
to \*inflate\* a scalar `ratm` object such that it can be, e.g.,
elementwise added. It is used after
[`upgrade_objects`](https://bfunovits.github.io/rationalmatrices/reference/upgrade_objects.md)
is called and therefore only relevant for

- [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),

- [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),

- [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),

- [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),

- [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)

In other words, one cannot inflate
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
or
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
objects.

## Usage

``` r
inflate_object(e, m, n)
```

## Arguments

- e:

  Object of type
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)

- m:

  Integer. Output dimension to which the scalar should be inflated.

- n:

  Integer. Output dimension to which the scalar should be inflated.

## Value

Object of same type but with \*inflated\* dimension

## Examples

``` r
(lp = test_lpolm(degree_max = 1, degree_min = -1))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>       z^-1 [,1]   z^0 [,1]   z^1 [,1]
#> [1,] -0.1742999 -0.3163715 -0.4930019
inflate_object(lp, 2, 2)
#> ( 2 x 2 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>       z^-1 [,1]       [,2]   z^0 [,1]       [,2]   z^1 [,1]       [,2]
#> [1,] -0.1742999 -0.1742999 -0.3163715 -0.3163715 -0.4930019 -0.4930019
#> [2,] -0.1742999 -0.1742999 -0.3163715 -0.3163715 -0.4930019 -0.4930019

(ss = test_stsp(s = 1))
#> statespace realization [1,1] with s = 1 states
#>            s[1]       u[1]
#> s[1]  0.1907785 -0.1064133
#> x[1] -1.1037450  1.0000000
inflate_object(ss, 2, 2)
#> statespace realization [2,2] with s = 2 states
#>            s[1]       s[2]       u[1]       u[2]
#> s[1]  0.1907785  0.0000000 -0.1064133 -0.1064133
#> s[2]  0.0000000  0.1907785 -0.1064133 -0.1064133
#> x[1] -1.1037450  0.0000000  1.0000000  1.0000000
#> x[2]  0.0000000 -1.1037450  1.0000000  1.0000000
```
