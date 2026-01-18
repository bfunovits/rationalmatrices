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
#>      z^0 [,1]  [,2]  [,3]    z^1 [,1]       [,2]       [,3]
#> [1,]        1     0     0  1.15842878 -0.6129769 -0.7270701
#> [2,]        0     1     0 -0.02718106  0.3336518  0.3520241
#> [3,]        0     0     1 -0.02928171  0.4855265 -1.7325593
#> right factor b(z):
#>        z^0 [,1]       [,2]   z^1 [,1]       [,2]
#> [1,]  0.4517525 -0.3768096  0.7161331 -1.1965556
#> [2,] -1.6380585  0.2512775 -0.6930620 -0.4940730
#> [3,]  1.2727707  2.0117546 -3.4547132  0.2029626
t(x)
#> ( 2 x 3 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>        z^0 [,1]       [,2]     [,3]   z^1 [,1]      [,2]       [,3]
#> [1,]  0.4517525 -1.6380585 1.272771  0.7161331 -0.693062 -3.4547132
#> [2,] -0.3768096  0.2512775 2.011755 -1.1965556 -0.494073  0.2029626
#> right factor c(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]        [,2]        [,3]
#> [1,]        1     0     0  1.1584288 -0.02718106 -0.02928171
#> [2,]        0     1     0 -0.6129769  0.33365185  0.48552650
#> [3,]        0     0     1 -0.7270701  0.35202409 -1.73255925
all.equal(x, t(t(x)))
#> [1] TRUE
```
