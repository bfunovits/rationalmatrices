# Constructor for Left Matrix Fraction Descriptions (LMFDs)

A Left Matrix Fraction Description (LMFD) of a rational matrix, \\x(z)\\
say, is a pair \\(a(z),b(z))\\ of polynomial matrices, such that \\x(z)
= a^{-1}(z) b(z)\\. The polynomial matrix \\a(z)\\ must be square and
invertible.

## Usage

``` r
lmfd(a, b)
```

## Arguments

- a, b:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects, or objects which may be coerced to a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object, via `x = polm(x)`. Either of the two arguments may be omitted.

## Value

An object of class `lmfd`.

## Details

Suppose that \\x(z)=a^{-1}(z) b(z)\\ is an \\(m,n)\\-dimensional matrix
and that \\a(z)\\ and \\b(z)\\ have degrees \\p\\ and \\q\\
respectively. The corresponding `lmfd` object stores the coefficients of
the polynomials \\a(z), b(z)\\ in an \\(m,m(p+1)+n(q+1))\\ dimensional
(real or complex valued) matrix together with an attribute
`order = c(m,n,p,q)` and a class attribute `c("lmfd", "ratm")`.

For a valid LMFD we require \\m\>0\\ and \\p\geq 0\\.

## See also

Useful methods and functions for the `lmfd` class are:

- [`test_lmfd`](https://bfunovits.github.io/rationalmatrices/reference/test_lmfd.md)
  generates random rational matrices in LMFD form.

- checks:
  [`is.lmfd`](https://bfunovits.github.io/rationalmatrices/reference/is.md),
  [`is.miniphase`](https://bfunovits.github.io/rationalmatrices/reference/check.md),
  [`is.stable`](https://bfunovits.github.io/rationalmatrices/reference/check.md)
  and
  [`is.coprime`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md).

- generic S3 methods: [`dim`](https://rdrr.io/r/base/dim.html),
  [`str`](https://rdrr.io/r/utils/str.html) and
  [`print`](https://rdrr.io/r/base/print.html).

- arithmetics:
  [`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md).

- matrix operations:
  [`bind`](https://bfunovits.github.io/rationalmatrices/reference/bind.md).

- extract the factors \\a(z)\\ and \\b(z)\\ with
  [`$.lmfd`](https://bfunovits.github.io/rationalmatrices/reference/extract.md).

- [`zeroes`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md),
  ...

## Examples

``` r
### (1 x 1) rational matrix x(z) = (1+z+z^2)^(-1) (3+2z+z^2)
lmfd(c(1,1,1), c(3,2,1)) %>% print(format = 'c')
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 2)
#> left factor a(z):
#>              [,1]
#> [1,]  1 + z + z^2 
#> right factor b(z):
#>               [,1]
#> [1,]  3 + 2z + z^2 

### (1 x 1) rational matrix x(z) = (3+2z+z^2)
lmfd(b = c(3,2,1)) %>% print(format = 'c')
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 0, q = 2)
#> left factor a(z):
#>       [,1]
#> [1,]     1 
#> right factor b(z):
#>               [,1]
#> [1,]  3 + 2z + z^2 

### (1 x 1) rational matrix x(z) = (1+z+z^2)^(-1)
lmfd(c(1,1,1)) %>% print(format = 'c')
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 0)
#> left factor a(z):
#>              [,1]
#> [1,]  1 + z + z^2 
#> right factor b(z):
#>       [,1]
#> [1,]     1 

### (2 x 3) rational matrix with degrees p=1, q=1
x = lmfd(array(rnorm(2*2*2), dim = c(2,2,2)), 
         array(rnorm(2*3*2), dim = c(2,3,2)))
is.lmfd(x)
#> [1] TRUE
dim(x)
#> m n p q 
#> 2 3 1 1 
str(x)
#> ( 2 x 3 ) left matrix fraction description with degrees (p = 1, q = 1)
print(x, digits = 2)
#> ( 2 x 3 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]    -0.75 -1.66     0.35 -0.11
#> [2,]    -0.73 -0.69     0.48 -0.99
#> right factor b(z):
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3]
#> [1,]    -0.20  0.96 -0.66     0.57  2.16  -0.2
#> [2,]     0.04  0.14 -1.08     2.08  1.10  -0.1

if (FALSE) { # \dontrun{
### the following calls to lmfd() throw an error 
lmfd() # no arguments!
lmfd(a = test_polm(dim = c(2,3), degree = 1))  # a(z) must be square 
lmfd(a = test_polm(dim = c(2,2), degree = -1)) # a(z) must have degree >= 0
} # }
```
