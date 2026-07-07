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
#>        z^0 [,1]       [,2]    z^1 [,1]       [,2]
#> [1,]  0.3404245 -1.5289587  0.74702859 -0.6395348
#> [2,] -0.4727025  0.2374253 -1.56251843 -0.8451957
#> [3,]  0.7087531 -1.3128142  0.07105336  0.6752447
#> right factor c(z):
#>      z^0 [,1]  [,2]   z^1 [,1]        [,2]   z^2 [,1]       [,2]
#> [1,]        1     0 -1.4522998 -0.33893587 0.04020439 -0.9984326
#> [2,]        0     1  0.9412061 -0.07557425 0.12430107  1.2333901
as.stsp(rr)
#> statespace realization [3,2] with s = 4 states
#>            s[1]        s[2]        s[3]       s[4]       u[1]       u[2]
#> s[1]  1.4522998  0.33893587 -0.04020439  0.9984326  1.0000000  0.0000000
#> s[2] -0.9412061  0.07557425 -0.12430107 -1.2333901  0.0000000  1.0000000
#> s[3]  1.0000000  0.00000000  0.00000000  0.0000000  0.0000000  0.0000000
#> s[4]  0.0000000  1.00000000  0.00000000  0.0000000  0.0000000  0.0000000
#> x[1]  2.6804923 -0.63970260  0.17636464  2.2256934  0.3404245 -1.5289587
#> x[2] -2.4724903 -0.98746832 -0.01050751 -0.7647996 -0.4727025  0.2374253
#> x[3]  2.3360041  0.81625159  0.13468922  2.3268542  0.7087531 -1.3128142
```
