# Construct Hankel Matrix from Impulse Response Coefficients

This (internal) helper function builds a Hankel matrix from a given
impulse response of dimension \\(m,n)\\ with \\l\\ lags. If the
parameter `Hsize=c(f,p)` is not given, then a default choice for the
number of block rows (\\f\\) and block columns (\\p\\) is made such that
\\f+p-1 = l\\, \\p\geq 1\\ and \\f \geq p+1\\.  
The function throws an error if the conditions \\p\geq 1\\, \\f \geq 2\\
and \\l \geq f+p-1\\ are not satisfied. In particular, this implies that
\\l \geq 2\\ must hold.

## Usage

``` r
pseries2hankel(obj, Hsize = NULL)
```

## Arguments

- obj:

  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  object or 3-D array with dimension \\(m,n,l+1)\\.

- Hsize:

  integer vector `c(f,p)`, number of block rows and block columns. If
  NULL a default choice is made.

## Value

Block Hankel matrix with dimension \\(fm, pn)\\ and attributes
`order=c(m,n,f,p)` and `k0`. The \\(m,n)\\-dimensional matrix `k0` is
the lag zero coefficient of the impulse response.
