# Check Properties of Rational Matrices

Check Properties of Rational Matrices

## Usage

``` r
is.stable(x, ...)

# S3 method for class 'ratm'
is.stable(x, ...)

is.miniphase(x, ...)

# S3 method for class 'ratm'
is.miniphase(x, ...)
```

## Arguments

- x:

  a rational matrix object, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  object.

- ...:

  not used.

## Value

Boolean. Note `is.stable`, `is.miniphase` return `NA` for objects of
type `pseries` and `zvalues`.

## Examples

``` r
a = test_polm(dim = c(2,2), deg = 1, random = TRUE)
b = test_polm(dim = c(2,3), deg = 1, random = TRUE)
c = lmfd(a,b)

is.stable(a)
#> [1] TRUE
is.stable(c)
#> [1] FALSE
is.stable(as.stsp(c))
#> [1] FALSE
is.stable(pseries(c))
#> [1] NA

is.miniphase(b[,1:2])
#> [1] FALSE
is.miniphase(as.stsp(c)[,1:2])
#> [1] FALSE

if (FALSE) { # \dontrun{
is.miniphase(b)
is.miniphase(c)
} # }
```
