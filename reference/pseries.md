# Power series Parameters

This function returns the coefficients of the power series expansion
(around \\z=0\\) of a rational matrix \$\$c(z) = c_0 + c_1 z + c_2 z^2 +
\cdots\$\$ Note that for a left matrix fraction description \\c(z) =
a(z)^{-1} b(z)\\ the matrix \\a(0)\\ must be an invertible matrix. If
the rational functions corresponds to a VARMA or statespace model, then
this sequence of coefficients is also called *impulse response
function*.

## Usage

``` r
pseries(obj, lag.max, ...)

# Default S3 method
pseries(obj, lag.max = 5, ...)

# S3 method for class 'polm'
pseries(obj, lag.max = 5, ...)

# S3 method for class 'lmfd'
pseries(obj, lag.max = 5, ...)

# S3 method for class 'rmfd'
pseries(obj, lag.max = 5, ...)

# S3 method for class 'stsp'
pseries(obj, lag.max = 5, ...)

# S3 method for class 'lpolm'
pseries(obj, ...)
```

## Arguments

- obj:

  (rational) matrix object, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object or an object which may be coerced to a polynomial matrix with
  `polm(obj)`. The default `S3` method first coerces the input argument
  `obj` to a `polm` object. If this fails an error is thrown.

- lag.max:

  (integer) maximum lag.

- ...:

  not used.

## Value

object of type `pseries`

## Details

A `pseries` object is simply a 3-dimensional (numeric or complex) array
of dimension c(m,n,lag.max+1) with a class attribute c('pseries',
'ratm'). The dimensions \\(m,n)\\ may also be zero.

## Examples

``` r
(obj = rmfd(c = polm(array(c(c(diag(2)), c(0.5, 0.25, 0.125, 0.5)), dim = c(2,2,2))), d = NULL))
#> ( 2 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 0
#> left factor d(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> right factor c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0     0.50 0.125
#> [2,]        0     1     0.25 0.500
(k = pseries(obj, lag.max = 4))
#> ( 2 x 2 ) impulse response with maximum lag = 4 
#>      lag=0 [,1]  [,2] lag=1 [,1]   [,2] lag=2 [,1]    [,2] lag=3 [,1]
#> [1,]          1     0      -0.50 -0.125    0.28125 0.12500 -0.1718750
#> [2,]          0     1      -0.25 -0.500    0.25000 0.28125 -0.1953125
#>             [,2] lag=4 [,1]      [,2]
#> [1,] -0.09765625  0.1103516 0.0703125
#> [2,] -0.17187500  0.1406250 0.1103516
pseries2rmfd(k)
#> $Xr
#> ( 2 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>      z^0 [,1]  [,2]     z^1 [,1]          [,2]
#> [1,]        1     0 0.000000e+00 -1.526557e-16
#> [2,]        0     1 5.551115e-17  0.000000e+00
#> right factor c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0     0.50 0.125
#> [2,]        0     1     0.25 0.500
#> 
#> $mu
#> [1] 1 1
#> 

(obj = lmfd(a = polm(array(c(c(diag(2)), c(0.5, 0.25, 0.125, 0.5)), dim = c(2,2,2)))))
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 0)
#> left factor a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0     0.50 0.125
#> [2,]        0     1     0.25 0.500
#> right factor b(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
(k = pseries(obj, lag.max = 4))
#> ( 2 x 2 ) impulse response with maximum lag = 4 
#>      lag=0 [,1]  [,2] lag=1 [,1]   [,2] lag=2 [,1]    [,2] lag=3 [,1]
#> [1,]          1     0      -0.50 -0.125    0.28125 0.12500 -0.1718750
#> [2,]          0     1      -0.25 -0.500    0.25000 0.28125 -0.1953125
#>             [,2] lag=4 [,1]      [,2]
#> [1,] -0.09765625  0.1103516 0.0703125
#> [2,] -0.17187500  0.1406250 0.1103516
pseries2lmfd(k)
#> $Xl
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0     0.50 0.125
#> [2,]        0     1     0.25 0.500
#> right factor b(z):
#>      z^0 [,1]  [,2]      z^1 [,1]         [,2]
#> [1,]        1     0 -1.110223e-16 2.775558e-17
#> [2,]        0     1  0.000000e+00 0.000000e+00
#> 
#> $nu
#> [1] 1 1
#> 
```
