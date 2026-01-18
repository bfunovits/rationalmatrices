# Derivative of a rational Matrix

Computes the derivative of a rational matrix \\k(z)\\ (with repect to
the complex variable \\z\\). Note that computing the derivative for an
impulse response object decreases the number of lags by one!

## Usage

``` r
derivative(obj, ...)

# S3 method for class 'lpolm'
derivative(obj, ...)

# S3 method for class 'polm'
derivative(obj, ...)

# S3 method for class 'lmfd'
derivative(obj, ...)

# S3 method for class 'rmfd'
derivative(obj, ...)

# S3 method for class 'stsp'
derivative(obj, ...)

# S3 method for class 'pseries'
derivative(obj, ...)

# S3 method for class 'zvalues'
derivative(obj, ...)
```

## Arguments

- obj:

  an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  or
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md).

- ...:

  not used

## Value

an object of the same class as the argument `obj`

## Examples

``` r
# create random (3 by 2) polynomial matrix with degree 2
K = test_polm(dim = c(3,2), degree = 2)
derivative(K)
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]      111   121      224   244
#> [2,]      211   221      424   444
#> [3,]      311   321      624   644
derivative(derivative(K))
#> ( 3 x 2 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]  [,2]
#> [1,]      224   244
#> [2,]      424   444
#> [3,]      624   644
derivative(derivative(derivative(K)))
#> ( 3 x 2 ) matrix polynomial with degree <= -1 

# note: computing the derivative of the impulse response 
# decreases "lag.max" by one!
all.equal(pseries(derivative(K)), derivative(pseries(K, lag.max = 6)))
#> [1] TRUE

# create statespace realization of a random (3 by 2) rational matrix
# with statespace dimension s = 4
K = test_stsp(dim = c(2,2), s = 4, bpoles = 1)
all.equal(pseries(derivative(K)), derivative(pseries(K, lag.max = 6)))
#> [1] TRUE

if (FALSE) { # \dontrun{
# 'lmfd' objects and 'zvalues' objects are not supported
derivative(test_lmfd(dim = c(3,3), degrees = c(1,1)))
derivative(zvalues(K))
} # }
```
