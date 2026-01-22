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
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]      [,2]       [,3]
#> [1,]        1     0     0 -0.7464070 0.3118272 -1.1678667
#> [2,]        0     1     0  0.8875486 1.4075345  0.8027295
#> [3,]        0     0     1 -0.2247582 1.3457063  0.1375144
#> right factor b(z):
#>        z^0 [,1]       [,2]  z^1 [,1]       [,2]
#> [1,] -0.1685889 -1.7743396 0.8989915  0.2191328
#> [2,] -1.6127394  0.6723864 0.6837201 -0.1998670
#> [3,]  1.3621800 -1.1277085 0.9352828  0.4877029
t(x)
#> ( 2 x 3 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>        z^0 [,1]       [,2]      [,3]  z^1 [,1]       [,2]      [,3]
#> [1,] -0.1685889 -1.6127394  1.362180 0.8989915  0.6837201 0.9352828
#> [2,] -1.7743396  0.6723864 -1.127709 0.2191328 -0.1998670 0.4877029
#> right factor c(z):
#>      z^0 [,1]  [,2]  [,3]   z^1 [,1]      [,2]       [,3]
#> [1,]        1     0     0 -0.7464070 0.8875486 -0.2247582
#> [2,]        0     1     0  0.3118272 1.4075345  1.3457063
#> [3,]        0     0     1 -1.1678667 0.8027295  0.1375144
all.equal(x, t(t(x)))
#> [1] TRUE
```
