# Polynomial Degree

Compute the polynomial degrees (of the elements) of a polynomial matrix.
Note that for a (scalar) polynomial with zero coefficients the degree is
set to \\-1\\.

## Usage

``` r
degree(x, which = c("elements", "rows", "columns", "matrix"))
```

## Arguments

- x:

  A polynomial matrix, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

- which:

  (character string) decides whether a matrix with the respectives
  degrees of the entries of the matrix, or a vector with the respective
  maximal degrees in each row or column, or simply the maximum degree of
  all elements of the polynomial matrix is computed.

## Value

The outcome depends on the parameter `which`:

- elements:

  A matrix with the degrees of the respective elements of the polynomial
  matrix.

- rows:

  A vector with the maximum degrees within each row.

- columns:

  A vector with the maximum degrees within each column.

- matrix:

  maximum of the degrees of the elements of the matrix.

## Details

The main advantage of setting the degree of a zero polynomial to \\-1\\
rather than `-Inf` is that indexing and assigning is programmatically
easier: E.g., `p[,,(deg+2):(max_deg+1)] = 0` (in obvious notation) also
works when `deg = -1`. When multiplying with zero polynomials, it does
not make a difference whether one needs to check for `deg = -Inf` or
`deg = -1`.

## Examples

``` r
x = polm(array(c(0,1,1,0,
                 0,0,1,0,
                 0,0,0,1,
                 0,0,0,0), dim = c(2,2,4)))
x
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]        0     1        0     1        0     0        0     0
#> [2,]        1     0        0     0        0     1        0     0
degree(x)
#>      [,1] [,2]
#> [1,]   -1    1
#> [2,]    0    2
degree(x, 'rows')
#> [1] 1 2
degree(x, 'columns')
#> [1] 0 2
degree(x, 'matrix')
#> [1] 2
```
