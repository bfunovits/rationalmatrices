# Compute left-Kronecker Indices for given Hankel Matrix

This (internal) helper functions determines the left-Kronecker indices
given the Hankel matrix of the impulse response coefficients. There are
no checks on parameters! In particular, note that the Hankel matrix must
have an attribute `order=c(m,n,f,p)` which describes the block size
\\(m,n)\\ and the number of block rows (\\f\\) and block columns
(\\p\\).

## Usage

``` r
hankel2nu(H, tol = sqrt(.Machine$double.eps))
```

## Arguments

- H:

  Block Hankel matrix, as computed e.g. by
  [`pseries2hankel`](https://bfunovits.github.io/rationalmatrices/reference/pseries2hankel.md).

- tol:

  tolerance parameter, used by [`qr`](https://rdrr.io/r/base/qr.html).

## Value

Integer vector with the Kronecker indices.
