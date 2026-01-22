# Construct an RMFD Representation from Impulse Response

This (helper) function constructs a right matrix fraction description of
a rational matrix from its impulse response. Of course, the impulse
response must contain sufficently many lags. The constructed RMFD is in
canonical form.

## Usage

``` r
pseries2rmfd(obj, Hsize = NULL, mu = NULL, tol = sqrt(.Machine$double.eps))
```

## Arguments

- obj:

  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  object or 3-D array with dimension \\(m,n,l+1)\\.

- Hsize:

  integer vector `c(f,p)`, number of block rows and block columns of the
  Hankel matrix which is used to construct the RMFD. If NULL a default
  choice is made.

- mu:

  integer vector with the right-Kronecker indices. If `NULL` then the
  right-Kronecker indices are computed via a QR decomposition of the
  transpose of the Hankel matrix, see
  [`qr`](https://rdrr.io/r/base/qr.html).

- tol:

  tolerance parameter, used by [`qr`](https://rdrr.io/r/base/qr.html).

## Value

List with components

- Xr:

  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
  object which contains the RMFD resresentation of the rational matrix
  (in echelon form).

- mu:

  integer vector with the right-Kronecker indices.

## Details

There are a number of restrictions on the dimension \\(m,n)\\ and number
of lags \\l\\ of the impulse response, the number of block rows (\\f\\),
block columns (\\p\\) of the Hankel matrix and the right-Kronecker
indices \\\mu_i\\:

We require that:

- \\m\>0\\, \\p\>0\\, \\f\>1\\,

- \\l \geq f+p-1\\ and

- \\\nu_i \<f\\.

If these restrictions are not satisfied an error is thrown.

## Examples

``` r
# generate a random RMFD object
m = 3
n = 2
p = 1
q = 1
cc = test_polm(dim = c(n,n), degree = p, random = TRUE)
dd = test_polm(dim = c(m,n), degree = q, random = TRUE)
Xr = rmfd(cc,dd)

# compute impulse response of this matrix
Xi = pseries(Xr, lag.max = 2*(max(p,q)+1))

# reconstruct a matrix from this impulse response
out = pseries2rmfd(Xi)
#> Warning: 'X' matrix was collinear
#> Warning: 'X' matrix was collinear
out
#> $Xr
#> ( 3 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>       z^0 [,1]       [,2]   z^1 [,1]       [,2]
#> [1,]  49.02065  39.684749  -14.82652  -12.10711
#> [2,] 148.21825 120.537069 -248.87162 -201.31285
#> [3,] -10.25606  -8.293895   41.56984   33.61174
#> right factor c(z):
#>      z^0 [,1]  [,2] z^1 [,1]     [,2]
#> [1,]        1     0 88.46234 70.67135
#> [2,]        0     1  0.00000  0.00000
#> 
#> $mu
#> [1] 1 1
#> 

# check that the lmfd object is a realization of the given impulse response
Xi1 = pseries(out$Xr, lag.max = 2*(max(p,q)+1))
all.equal(Xi, Xi1)
#> [1] "Mean relative difference: 0.0002085268"
```
