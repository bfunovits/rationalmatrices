# Tools related to Kronecker indices

The Kronecker indices \\(\nu_1,\ldots,\nu_m)\\ describe a (nice) basis
of the row space of the Hankel matrix of the impulse response
coefficients. These indices are e.g. used to construct a (unique) LMFD
representation in echelon canonical form for a given impulse response.

## Usage

``` r
basis2nu(basis, m)

nu2basis(nu)

pseries2nu(obj, Hsize = NULL, tol = sqrt(.Machine$double.eps))
```

## Arguments

- basis:

  s-dimensional (integer) vector which contains the indices of the basis
  rows of the Hankel matrix.

- m:

  (integer) the number of rows of the underlying rational matrix.

- nu:

  vector with the Kronecker indices, i.e. an m-dimensional integer
  vector.

- obj:

  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  object or 3-D array with dimension \\(m,n,l+1)\\. This object
  represents the impulse response of a rational matrix.

- Hsize:

  integer vector `c(f,p)`, number of block rows and block columns of the
  Hankel matrix of the impulse response coefficients. If NULL a default
  choice is made.

- tol:

  tolerance parameter, used by [`qr`](https://rdrr.io/r/base/qr.html).

## Value

- basis2nu:

  returns the Kronecker indices (`nu`) for given indices of the basis
  rows.

- nu2basis:

  returns the indices of the basis rows for given Kronecker indices
  `nu`.

- pseries2nu:

  determines the Kronecker indices for given impulse response
  coefficients. The function uses a QR decomposition (with pivoting) of
  the transpose of the Hankel matrix to determine the rank and the basis
  for its row space. See [`qr`](https://rdrr.io/r/base/qr.html).

## Details

The function `pseries2nu` first constructs a Hankel matrix of the
impulse response coefficients with \\f\\ and \\p\\ block rows and
columns respectively. Then a (nice) basis for the row space is computed
via a QR decomposition (with pivoting) of the transposed Hankel matrix.
See [`qr`](https://rdrr.io/r/base/qr.html). If the size of the Hankel
matrix is not specified (the parameter `Hsize=c(f,p)` is missing), then
a default choice is made such that \\f+p-1 = l\\, \\p\geq 1\\ and \\f
\geq p+1\\ holds.  
The function `pseries2nu` throws an error if the conditions \\p\geq 1\\,
\\f \geq 2\\ and \\l \geq f+p-1\\ are not satisfied. In particular, this
implies that \\l \geq 2\\ must hold.

## References

Hannan EJ, Deistler M (2012). *The Statistical Theory of Linear
Systems*, Classics in Applied Mathematics. SIAM, Philadelphia.
Originally published: John Wiley & Sons, New York, 1988.

## Examples

``` r
basis = c(1,2,3,4,6,7)     # suppose rows 1,2,3,4,6,7 of the Hankel matrix form a basis
nu = basis2nu(basis, m=3)  # compute the corresponding Kronecker index
print(nu)
#> [1] 3 1 2
all.equal(basis,nu2basis(nu))  # nu2basis(nu) returns the indices of the basis rows (basis)
#> [1] TRUE

# generate random rational matrix in statspace form
# make sure that the matrix is stable
m = 3
n = 2
s = 7
A = matrix(rnorm(s*s), nrow = s, ncol = s)
A = A / (1.1 * max(abs(eigen(A, only.values = TRUE)$values)))
x = stsp(A, B = matrix(rnorm(s*n), nrow = s, ncol = n),
            C = matrix(rnorm(s*m), nrow = m, ncol = s),
            D = diag(1, nrow = m, ncol = n))
k = pseries(x, lag.max = 20)

# compute the Kronecker indices of this  rational matrix
pseries2nu(k)
#> [1] 3 2 2

if (FALSE) { # \dontrun{
# Suppose the rational matrix has dimension m=2. Then the rows 1,2,5
# do not form a "nice" basis for the row space of the Hankel matrix.
# Therefore "basis2nu" stops with an error message.
basis2nu(c(1,2,5), m=2)
} # }
```
