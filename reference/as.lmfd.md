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

# S3 method for class 'stsp'
as.lmfd(
  obj,
  method = c("echelon"),
  lag.max = NULL,
  tol = sqrt(.Machine$double.eps),
  ...
)
```

## Arguments

- obj:

  object

- method:

  character string

- ...:

  optional additional arguments

- lag.max:

  Integer. Number of lags for the impulse response computation. Defaults
  to `max(2 * s, 10)` where `s` is the state dimension of `obj`. Must be
  large enough for the Hankel matrix to resolve the Kronecker indices.

- tol:

  Tolerance for the rank decision in the Kronecker-index computation
  (passed to [`qr`](https://rdrr.io/r/base/qr.html)). Default:
  `sqrt(.Machine$double.eps)`.

## Value

object of class
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
