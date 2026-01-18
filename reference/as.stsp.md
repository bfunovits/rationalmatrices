# Coerce to Statespace Realization

The function `as.stsp.pseries` calls
[`pseries2stsp`](https://bfunovits.github.io/rationalmatrices/reference/pseries2stsp.md)
with default parameters. Of course the
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
object must contain sufficiently many lags.

## Usage

``` r
as.stsp(obj, ...)

# S3 method for class 'lpolm'
as.stsp(obj, ...)

# S3 method for class 'polm'
as.stsp(obj, ...)

# S3 method for class 'lmfd'
as.stsp(obj, ...)

# S3 method for class 'rmfd'
as.stsp(obj, ...)

# S3 method for class 'pseries'
as.stsp(obj, method = c("balanced", "echelon"), ...)
```

## Arguments

- obj:

  object

- ...:

  optional additional parameters

- method:

  character string

## Value

object of class
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md).

## Examples

``` r
(rr = test_rmfd(dim = c(3,2), degrees = c(2,1)))
#> ( 3 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 2, deg(d(z)) = q = 1
#> left factor d(z):
#>         z^0 [,1]       [,2]   z^1 [,1]      [,2]
#> [1,] -1.56251843 -0.8451957 -1.6865047  1.100190
#> [2,]  0.07105336  0.6752447 -0.9028149  1.203768
#> [3,] -0.63953477  1.1533758  1.3176337 -1.431271
#> right factor c(z):
#>      z^0 [,1]  [,2]  z^1 [,1]       [,2]   z^2 [,1]       [,2]
#> [1,]        1     0 1.2333901 -0.4727025 -1.5289587 -1.3128142
#> [2,]        0     1 0.3404245  0.7087531  0.2374253  0.7470286
as.stsp(rr)
#> statespace realization [3,2] with s = 4 states
#>            s[1]       s[2]        s[3]       s[4]        u[1]       u[2]
#> s[1] -1.2333901  0.4727025  1.52895871  1.3128142  1.00000000  0.0000000
#> s[2] -0.3404245 -0.7087531 -0.23742535 -0.7470286  0.00000000  1.0000000
#> s[3]  1.0000000  0.0000000  0.00000000  0.0000000  0.00000000  0.0000000
#> s[4]  0.0000000  1.0000000  0.00000000  0.0000000  0.00000000  0.0000000
#> x[1]  0.5284153  0.9606185 -2.18835529 -1.4199111 -1.56251843 -0.8451957
#> x[2] -1.2203213  0.7587732 -0.05168255 -0.4111472  0.07105336  0.6752447
#> x[3]  1.7137922 -2.5510391 -1.25166291 -1.7011950 -0.63953477  1.1533758
```
