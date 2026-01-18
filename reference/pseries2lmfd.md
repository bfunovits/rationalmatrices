# Construct a LMFD Representation from Impulse Response

This (helper) function constructs a left matrix fraction description of
a rational matrix from its impulse response. Of course the impulse
response must contain sufficiently many lags. The constructed LMFD is in
*echelon canonical* form.

## Usage

``` r
pseries2lmfd(obj, Hsize = NULL, nu = NULL, tol = sqrt(.Machine$double.eps))
```

## Arguments

- obj:

  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  object or 3-D array with dimension \\(m,n,l+1)\\.

- Hsize:

  integer vector `c(f,p)`, number of block rows and block columns of the
  Hankel matrix which is used to construct the LMFD. If NULL a default
  choice is made.

- nu:

  integer vector with the Kronecker indices. If `NULL` then the
  Kronecker indices are computed via a QR decomposition of the transpose
  of the Hankel matrix, see [`qr`](https://rdrr.io/r/base/qr.html).

- tol:

  tolerance parameter, used by [`qr`](https://rdrr.io/r/base/qr.html).

## Value

List with components

- Xl:

  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  object which contains the LMFD representation of the rational matrix
  (in echelon form).

- nu:

  integer vector with the Kronecker indices.

## Details

There are a number of restrictions on the dimension \\(m,n)\\ and number
of lags \\l\\ of the impulse response, the number of block rows (\\f\\),
block columns (\\p\\) of the Hankel matrix and the Kronecker indices
\\\nu_i\\:

We require that: \\m\>0\\, \\p\>0\\, \\f\>1\\, \\l \geq f+p-1\\ and
\\\nu_i \<f\\. If these restrictions are not satisfied an error is
thrown.

## Examples

``` r
# generate a random LMFD object
m = 3
n = 2
p = 1
q = 1
a = test_polm(dim = c(m,m), degree = p, random = TRUE, bzeroes = 1)
b = test_polm(dim = c(m,n), degree = q, random = TRUE, bzeroes = 1)
Xl = lmfd(a,b)

# compute impulse response of this matrix
Xi = pseries(Xl, lag.max = 2*(max(p,q)+1))

# reconstruct a matrix from this impulse response
out = pseries2lmfd(Xi)
out
#> $Xl
#> ( 3 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2]  [,3]     z^1 [,1]       [,2]       [,3]
#> [1,]        1     0     0 -0.006911964  0.7127661 -1.5267106
#> [2,]        0     1     0 -0.003477361  0.2490580 -0.4782623
#> [3,]        0     0     1  0.049013059 -0.6321853 -0.3212550
#> right factor b(z):
#>       z^0 [,1]        [,2]   z^1 [,1]       [,2]
#> [1,] 1.8271897 -0.87821457 -1.1016513 0.57336537
#> [2,] 0.7698252 -0.08715286 -0.1636610 0.15127154
#> [3,] 0.2278321  1.03854642 -0.1061596 0.04172095
#> 
#> $nu
#> [1] 1 1 1
#> 

# check that the lmfd object is a realization of the given impulse response
Xi1 = pseries(out$Xl, lag.max = 2*(max(p,q)+1))
all.equal(Xi, Xi1)
#> [1] TRUE
```
