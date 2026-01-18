# Constructor for Right Matrix Fraction Descriptions (RMFDs)

A Right Matrix Fraction Description (RMFD) of a rational matrix,
\\k(z)\\ say, is pair \\(c(z),d(z))\\ of polynomial matrices, such that
\\k(z) = d(z) c^{-1}(z)\\. The polynomial matrix \\c(z)\\ must be square
and invertible.

## Usage

``` r
rmfd(c = NULL, d = NULL)
```

## Arguments

- c, d:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects, or objects which may be coerced to a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object, via `x = polm(x)`. Either of the two arguments may be omitted.

## Value

An object of class `rmfd`.

## Details

Suppose that \\k(z) = d(z) c^{-1}(z)\\ is an \\(m,n)\\-dimensional
matrix and that \\c(z)\\ and \\d(z)\\ have degrees \\p\\ and \\q\\
respectively. The corresponding `rmfd` object stores the coefficients of
the polynomials \\c(z), d(z)\\ in an \\(n(p+1)+m(q+1), n)\\ dimensional
(real or complex valued) matrix. The `rmfd` object also stores the
attribute `order = c(m,n,p,q)` and a class attribute
`c("lmfd", "ratm")`.

For a valid RMFD we require \\m\>0\\ and \\p\geq 0\\.

## See also

Useful methods and functions for the `rmfd` class are:

- [`test_rmfd`](https://bfunovits.github.io/rationalmatrices/reference/test_rmfd.md)
  generates random rational matrices in RMFD form.

- checks:
  [`is.rmfd`](https://bfunovits.github.io/rationalmatrices/reference/is.md),
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

- extract the factors \\c(z)\\ and \\d(z)\\ with
  [`$.rmfd`](https://bfunovits.github.io/rationalmatrices/reference/extract.md).

- [`zeroes`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md),
  ...

## Examples

``` r
### (1 x 1) rational matrix k(z) =  (3+2z+z^2) (1+z+z^2)^(-1)
rmfd(c(1,1,1), c(3,2,1)) %>% print(format = 'c')
#> ( 1 x 1 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 2, deg(d(z)) = q = 2
#> left factor d(z):
#>               [,1]
#> [1,]  3 + 2z + z^2 
#> right factor c(z):
#>              [,1]
#> [1,]  1 + z + z^2 

### (1 x 1) rational matrix k(z) = (3+2z+z^2)
rmfd(c = c(3,2,1)) %>% print(format = 'c')
#> ( 1 x 1 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 2, deg(d(z)) = q = 0
#> left factor d(z):
#>       [,1]
#> [1,]     1 
#> right factor c(z):
#>               [,1]
#> [1,]  3 + 2z + z^2 

### (1 x 1) rational matrix k(z) = (1+z+z^2)^(-1)
rmfd(c(1,1,1)) %>% print(format = 'c')
#> ( 1 x 1 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 2, deg(d(z)) = q = 0
#> left factor d(z):
#>       [,1]
#> [1,]     1 
#> right factor c(z):
#>              [,1]
#> [1,]  1 + z + z^2 

### (3 x 2) rational matrix with degrees p=1, q=1
x = rmfd(array(rnorm(2*2*2), dim = c(2,2,2)), 
         array(rnorm(3*2*2), dim = c(3,2,2)))
is.rmfd(x)
#> [1] TRUE
dim(x)
#> m n p q 
#> 3 2 1 1 
str(x)
#> ( 3 x 2 ) right matrix fraction description with degrees (deg(c(z)) = p = 1, deg(d(z)) = q = 1)
print(x, digits = 2)
#> ( 3 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]     0.05  0.52    -1.29  2.07
#> [2,]     1.07 -0.20     1.07  0.46
#> [3,]     1.22 -1.61     1.83 -1.13
#> right factor c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]     1.52 -0.99     2.50 -0.80
#> [2,]    -0.18 -0.18     1.37  0.77

xr = test_rmfd(dim = c(3,2), degree = c(2,2))

if (FALSE) { # \dontrun{
### the following calls to rmfd() throw an error 
rmfd() # no arguments!
rmfd(c = test_polm(dim = c(2,3), degree = 1))  # c(z) must be square 
rmfd(c = test_polm(dim = c(2,2), degree = -1)) # c(z) must have degree >= 0
} # }
```
