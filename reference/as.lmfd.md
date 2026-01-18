# Coerce to Left Matrix Fraction Description

The function `as.lmfd.pseries` calls
[`pseries2lmfd`](https://bfunovits.github.io/rationalmatrices/reference/pseries2lmfd.md)
with default parameters. Of course the
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
object must contain sufficiently many lags.

## Usage

``` r
as.lmfd(obj, method, ...)

# S3 method for class 'pseries'
as.lmfd(obj, method, ...)
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
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
