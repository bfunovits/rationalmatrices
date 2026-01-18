# Frequency Response Function

This function evaluates a rational matrix at one given (complex)
arguments. It is the non-plural version of
[`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md).

## Usage

``` r
zvalue(obj, z, ...)

# S3 method for class 'lpolm'
zvalue(obj, z = NULL, ...)

# S3 method for class 'polm'
zvalue(obj, z = NULL, ...)

# S3 method for class 'lmfd'
zvalue(obj, z = NULL, ...)

# S3 method for class 'rmfd'
zvalue(obj, z = NULL, ...)

# S3 method for class 'stsp'
zvalue(obj, z = NULL, ...)
```

## Arguments

- obj:

  (rational) matrix object, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object or an object which may be coerced to a polynomial matrix with
  `polm(obj)`. The default `S3` method first coerces the input argument
  `obj` to a `polm` object. If this fails an error is thrown.

- z:

  (numeric or complex) vector of length one at which to evaluate the
  rational matrix.

- ...:

  optional additional parameters

## Value

Matrix of dimensions of input object
