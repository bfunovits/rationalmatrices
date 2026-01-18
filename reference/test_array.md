# Create Test Array

The helper function `test_array` creates arrays of given dimension.

## Usage

``` r
test_array(dim, random = FALSE, dimnames = FALSE)
```

## Arguments

- dim:

  Integer vector.

- random:

  Boolean. Default set to FALSE. If TRUE, the elements are randomly
  generated (standard normal).

- dimnames:

  (logical) decides whether `dimnames` should be created.

## Value

Real valued array of dimension `dim`.

## Details

The function `test_array` creates an array with either

- entries of the form `x[i,j,...] = i*10^(n-1) + j*10^(n-2) + ...`,
  where `n` is the number of dimensions,

- or randomly generated (standard normal) entries.

This function is mainly used to test operations on arrays (like
[`btoeplitz`](https://bfunovits.github.io/rationalmatrices/reference/btoeplitz.md),
[`bhankel`](https://bfunovits.github.io/rationalmatrices/reference/bhankel.md),
[`bmatrix`](https://bfunovits.github.io/rationalmatrices/reference/bmatrix.md)
and
[`dbind`](https://bfunovits.github.io/rationalmatrices/reference/dbind.md)).

## Examples

``` r
test_array(dim = 5)
#> [1] 1 2 3 4 5
test_array(dim = c(2,4), dimnames = TRUE)
#>      B
#> A     B=1 B=2 B=3 B=4
#>   A=1  11  12  13  14
#>   A=2  21  22  23  24
test_array(dim = c(2,4,3), dimnames = FALSE)
#> , , 1
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]  111  121  131  141
#> [2,]  211  221  231  241
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]  112  122  132  142
#> [2,]  212  222  232  242
#> 
#> , , 3
#> 
#>      [,1] [,2] [,3] [,4]
#> [1,]  113  123  133  143
#> [2,]  213  223  233  243
#> 

if (FALSE) { # \dontrun{
# the examples below throw an error
test_array(dim = c())
test_array(dim = c(2,-1,2))
} # }
```
