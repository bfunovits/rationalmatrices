# Block matrices

This helper function coerces an array to a (block) matrix.

## Usage

``` r
bmatrix(x, rows = NULL, cols = NULL)
```

## Arguments

- x:

  Vector, matrix or array. Vectors are coerced to one dimensional arrays
  and matrices are of course treated as 2-dimensional arrays.

- rows, cols:

  integer vectors. These two vectors define a partition of the
  "dimensions" `(1,...,n)`, where `n` is the number of dimensions of `x`
  (i.e. `length(dim(x))`). If either of the two is missing, then the
  complement is used. At least one of the arguments "rows" and "cols"
  has to be given.

## Value

matrix

## Examples

``` r
x = 1:4
bmatrix(x, rows = 1, cols = integer(0)) # returns an (4,1) matrix
#>      [,1]
#> [1,]    1
#> [2,]    2
#> [3,]    3
#> [4,]    4
bmatrix(x, cols = 1, rows = NULL)       # returns an (1,4) matrix
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    2    3    4

x = test_array(dim = c(2,3))
bmatrix(x, cols = 2)    # returns x    (is equivalent to bmatrix(x, rows = 1, cols = 2))
#>      [,1] [,2] [,3]
#> [1,]   11   12   13
#> [2,]   21   22   23
bmatrix(x, rows = 2)    # returns t(x) (is equivalent to bmatrix(x, rows = 2, cols = 1))
#>      [,1] [,2]
#> [1,]   11   21
#> [2,]   12   22
#> [3,]   13   23
bmatrix(x, rows = integer(0))   # returns an (1,6) matrix
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]   11   21   12   22   13   23

x = test_array(dim = c(2,3,4))
bmatrix(x, rows = 1)   # is equivalent to: bmatrix(x, rows = 1, cols = c(2,3))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#> [1,]  111  121  131  112  122  132  113  123  133   114   124   134
#> [2,]  211  221  231  212  222  232  213  223  233   214   224   234
bmatrix(x, rows = 1, cols = c(3,2))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#> [1,]  111  112  113  114  121  122  123  124  131   132   133   134
#> [2,]  211  212  213  214  221  222  223  224  231   232   233   234
bmatrix(x, cols = 2)   # is equivalent to: bmatrix(x, rows = c(1,3), cols = 2)
#>      [,1] [,2] [,3]
#> [1,]  111  121  131
#> [2,]  211  221  231
#> [3,]  112  122  132
#> [4,]  212  222  232
#> [5,]  113  123  133
#> [6,]  213  223  233
#> [7,]  114  124  134
#> [8,]  214  224  234
bmatrix(x, rows = 1:3) # is equivalent to: bmatrix(x, cols = integer(0))
#>       [,1]
#>  [1,]  111
#>  [2,]  211
#>  [3,]  121
#>  [4,]  221
#>  [5,]  131
#>  [6,]  231
#>  [7,]  112
#>  [8,]  212
#>  [9,]  122
#> [10,]  222
#> [11,]  132
#> [12,]  232
#> [13,]  113
#> [14,]  213
#> [15,]  123
#> [16,]  223
#> [17,]  133
#> [18,]  233
#> [19,]  114
#> [20,]  214
#> [21,]  124
#> [22,]  224
#> [23,]  134
#> [24,]  234
bmatrix(x, rows = c(3,1,2))
#>       [,1]
#>  [1,]  111
#>  [2,]  112
#>  [3,]  113
#>  [4,]  114
#>  [5,]  211
#>  [6,]  212
#>  [7,]  213
#>  [8,]  214
#>  [9,]  121
#> [10,]  122
#> [11,]  123
#> [12,]  124
#> [13,]  221
#> [14,]  222
#> [15,]  223
#> [16,]  224
#> [17,]  131
#> [18,]  132
#> [19,]  133
#> [20,]  134
#> [21,]  231
#> [22,]  232
#> [23,]  233
#> [24,]  234

if (FALSE) { # \dontrun{
# the examples below throw an error
bmatrix(x, rows = 1, cols = 2)
bmatrix(x, rows = c(1,2), cols = c(2,3))
bmatrix(x, rows = c(1,2,1), cols = 3)
} # }
```
