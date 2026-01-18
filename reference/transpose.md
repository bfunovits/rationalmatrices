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
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]       [,2]      [,3]
#> [1,]        1     0     0  0.4881914  0.3722990 0.2515316
#> [2,]        0     1     0 -1.3334324 -1.3931584 1.4003740
#> [3,]        0     0     1 -0.6412660  0.4756471 1.2336156
#> right factor b(z):
#>        z^0 [,1]       [,2]    z^1 [,1]       [,2]
#> [1,]  0.8713687  1.5276530 -0.40423152  0.4698279
#> [2,] -0.7459840 -0.1158758 -0.01323844 -1.9713205
#> [3,]  0.2781470  1.4555469 -0.59943140  0.3595583
t(x)
#> ( 2 x 3 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>       z^0 [,1]       [,2]     [,3]   z^1 [,1]        [,2]       [,3]
#> [1,] 0.8713687 -0.7459840 0.278147 -0.4042315 -0.01323844 -0.5994314
#> [2,] 1.5276530 -0.1158758 1.455547  0.4698279 -1.97132047  0.3595583
#> right factor c(z):
#>      z^0 [,1]  [,2]  [,3]  z^1 [,1]      [,2]       [,3]
#> [1,]        1     0     0 0.4881914 -1.333432 -0.6412660
#> [2,]        0     1     0 0.3722990 -1.393158  0.4756471
#> [3,]        0     0     1 0.2515316  1.400374  1.2336156
all.equal(x, t(t(x)))
#> [1] TRUE
```
