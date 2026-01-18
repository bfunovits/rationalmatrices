# Transform between Linear Index and Matrix Indices

These functions are as in MATLAB. `ind2sub()` transforms a linear index
to a row and column index for a matrix of given size. `sub2ind()`
transforms a matrix index, (row, col) to a a linear index (in terms of
columns).

## Usage

``` r
ind2sub(dim, ind)

sub2ind(dim, row, col)
```

## Arguments

- dim:

  Integer vector of size 2. Matrix dimensions.

- ind:

  Integer. Linear index.

- row:

  Integer. Row index.

- col:

  Integer. Column index.

## Value

`ind2sub()` returns for given linear index, a matrix index (row, col).
`sub2ind()` returns for a given matrix index (row, col), a linear index
(column-major).

## Examples

``` r
A = matrix(1:(3*4), 3, 4)
A
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    4    7   10
#> [2,]    2    5    8   11
#> [3,]    3    6    9   12

ind2sub(c(3,4), 7)
#> [1] 1 3
sub2ind(c(3,4), 2, 3)
#> [1] 8
```
