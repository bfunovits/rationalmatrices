# Rational Matrix Transpose

Compute the transpose of a rational matrix \\x(z)\\.

## Usage

``` r
# S3 method for class 'polm'
t(x)

# S3 method for class 'lpolm'
t(x)

# S3 method for class 'lmfd'
t(x)

# S3 method for class 'rmfd'
t(x)

# S3 method for class 'stsp'
t(x)

# S3 method for class 'pseries'
t(x)

# S3 method for class 'zvalues'
t(x)
```

## Arguments

- x:

  rational matrix object, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  object.

## Value

A rational matrix object, which represents the transposed rational
matrix \\x'(z)\\. The output is of the same class as the input `x`
unless `x` is an
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
or an
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
object: The transposition of an `lmfd` object is an `rmfd` object, and
vice versa.

## Examples

``` r
x = test_polm(dim = c(2,3), degree = 3, random = TRUE)
all.equal(pseries(t(x)), t(pseries(x)))
#> [1] TRUE

x = test_stsp(dim = c(3,2), s = 1)
all.equal(zvalues(t(x)), t(zvalues(x)))
#> [1] TRUE

# the transpose of an LMFD object is RMFD
(x = test_lmfd(dim = c(3,2), degrees = c(1,1)))
#> ( 3 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]       [,2]       [,3]
#> [1,]        1     0     0 0.04219385  1.2320383  2.1414703
#> [2,]        0     1     0 0.92030470  0.3614530 -1.0347436
#> [3,]        0     0     1 0.07766999 -0.2355909 -0.4225692
#> right factor b(z):
#>        z^0 [,1]       [,2]   z^1 [,1]        [,2]
#> [1,] -0.8148421  0.1449392 -0.8626400  0.02042165
#> [2,] -0.2761971 -0.2467375 -0.1293788 -0.24340699
#> [3,] -0.2846398  1.0116738 -1.2438879  0.38636828
t(x)
#> ( 2 x 3 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>        z^0 [,1]       [,2]       [,3]    z^1 [,1]       [,2]       [,3]
#> [1,] -0.8148421 -0.2761971 -0.2846398 -0.86263999 -0.1293788 -1.2438879
#> [2,]  0.1449392 -0.2467375  1.0116738  0.02042165 -0.2434070  0.3863683
#> right factor c(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]       [,2]        [,3]
#> [1,]        1     0     0 0.04219385  0.9203047  0.07766999
#> [2,]        0     1     0 1.23203830  0.3614530 -0.23559088
#> [3,]        0     0     1 2.14147034 -1.0347436 -0.42256919
all.equal(x, t(t(x)))
#> [1] TRUE
```
