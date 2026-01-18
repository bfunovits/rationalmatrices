# Helper for Extracting indices

Linear indices for extracting elements of a matrix are obtained,
depending on the dimension of the matrix, number of arguments `n_args`,
the missingness `is_missing` of these arguments, and their contents `i`
and `j`.

## Usage

``` r
extract_matrix_(m, n, n_args, is_missing, i, j)
```

## Arguments

- m:

  Integer. Output dimension.

- n:

  Integer. Input dimension.

- n_args:

  Integer. One, or two. Corresponds to the case where one set of indices
  or two sets are given. Corresponds to x\[i\] and x\[i,j\]

- is_missing:

  Boolean vector of dimension `n_args`, i.e. one or two. TRUE when this
  function is used fo x\[\] and x\[,\]

- i:

  First set of indices. Can be integers or booleans.

- j:

  Second set of indices. Like above.

## Value

Indices (integers) which will be further used to extract the relevant
elements

## Details

No input checks are performed because this is an internal function.
Consequently, some input combination are contradictory without alerting
the user.

## Examples

``` r
# x[,]
extract_matrix_(6, 5, 
                n_args = 2,
                is_missing = c(TRUE, TRUE),
                i = 1:3,
                j = c(2,5))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    7   13   19   25
#> [2,]    2    8   14   20   26
#> [3,]    3    9   15   21   27
#> [4,]    4   10   16   22   28
#> [5,]    5   11   17   23   29
#> [6,]    6   12   18   24   30
                
# x[i,j]
extract_matrix_(6, 5, 
                n_args = 2,
                is_missing = c(FALSE, FALSE),
                i = 1:3,
                j = c(2,5))
#>      [,1] [,2]
#> [1,]    7   25
#> [2,]    8   26
#> [3,]    9   27
                
# x[,j]
extract_matrix_(6, 5, 
                n_args = 2,
                is_missing = c(TRUE, FALSE),
                i = 1:3,
                j = c(2,5))
#>      [,1] [,2]
#> [1,]    7   25
#> [2,]    8   26
#> [3,]    9   27
#> [4,]   10   28
#> [5,]   11   29
#> [6,]   12   30
                
# x[i,]
extract_matrix_(6, 5, 
                n_args = 2,
                is_missing = c(FALSE, TRUE),
                i = 1:3,
                j = c(2,5))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    7   13   19   25
#> [2,]    2    8   14   20   26
#> [3,]    3    9   15   21   27
                
# x[i], j is ignored if available (doesn't happen because it will not be called in this way)
extract_matrix_(6, 5, 
                n_args = 1,
                is_missing = c(FALSE, TRUE),
                i = 1:10,
                j = c(2,5))
#>       [,1]
#>  [1,]    1
#>  [2,]    2
#>  [3,]    3
#>  [4,]    4
#>  [5,]    5
#>  [6,]    6
#>  [7,]    7
#>  [8,]    8
#>  [9,]    9
#> [10,]   10
                
# x[i], j is ignored if available (doesn't happen because it will not be called in this way)
extract_matrix_(6, 5, 
                n_args = 1,
                is_missing = c(FALSE, TRUE),
                i = c(3,6,7,1),
                j = c(2,5))
#>      [,1]
#> [1,]    3
#> [2,]    6
#> [3,]    7
#> [4,]    1
```
