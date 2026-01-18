# Coerce to Right Matrix Fraction Description

The function `as.rmfd.pseries` calls
[`pseries2rmfd`](https://bfunovits.github.io/rationalmatrices/reference/pseries2rmfd.md)
with default parameters. Of course the
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
object must contain sufficiently many lags.

## Usage

``` r
as.rmfd(obj, method, ...)

# S3 method for class 'pseries'
as.rmfd(obj, method, ...)
```

## Arguments

- obj:

  object

- method:

  character string

- ...:

  optional additional arguments

## Value

object of class
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
