# Block Hankel matrix

Construct a block Hankel matrix with (block) entries from a
3-dimensional array R. The \\(i,j)\\-th block of the Hankel matrix is
`R[,,i+j-1]`.

## Usage

``` r
bhankel(R, d = NULL)
```

## Arguments

- R:

  3-dimensional array, matrix or vector. A vector of length \\k\\ is
  coerced to a \\(1,1,k)\\-dimensional array and a \\(p,q)\\ matrix is
  treated as an array of dimension \\(p,q,1)\\.

- d:

  determines the number of block rows and columns. Suppose that `R` is
  an array of size \\(p,q,k)\\. If \\d=(m,n)\\ then `bhankel` returns a
  block Hankel matrix with \\m\\ block rows and \\n\\ block columns. If
  \\d=m\\ then a Hankel matrix with \\m\\ block rows and
  \\n=\max(k+1-m,1)\\ block columns is returned. In the default case
  `d = NULL` the number of block rows is \\m=(k+1)/2\\ for odd \\k\\ and
  \\m=(k/2+1)\\ else and the number of block columns is set to
  \\n=\max(k+1-m,1)\\. In the case \\(m+n-1)\>k\\ the array `R` is
  padded with zeros.

## Value

Matrix of size \\(pm\times qn)\\.

## Examples

``` r
bhankel(1:6)
#>      [,1] [,2] [,3]
#> [1,]    1    2    3
#> [2,]    2    3    4
#> [3,]    3    4    5
#> [4,]    4    5    6
bhankel(1:6, d = 3)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    2    3    4
#> [2,]    2    3    4    5
#> [3,]    3    4    5    6
bhankel(letters[1:6], d = c(3,7))  # note the "zero padding"
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
#> [1,] "a"  "b"  "c"  "d"  "e"  "f"  "0" 
#> [2,] "b"  "c"  "d"  "e"  "f"  "0"  "0" 
#> [3,] "c"  "d"  "e"  "f"  "0"  "0"  "0" 
bhankel(test_array(dim = c(2,2,6)))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  111  121  112  122  113  123
#> [2,]  211  221  212  222  213  223
#> [3,]  112  122  113  123  114  124
#> [4,]  212  222  213  223  214  224
#> [5,]  113  123  114  124  115  125
#> [6,]  213  223  214  224  215  225
#> [7,]  114  124  115  125  116  126
#> [8,]  214  224  215  225  216  226
bhankel(test_array(dim = c(3,2,6)), d = 3)
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#>  [1,]  111  121  112  122  113  123  114  124
#>  [2,]  211  221  212  222  213  223  214  224
#>  [3,]  311  321  312  322  313  323  314  324
#>  [4,]  112  122  113  123  114  124  115  125
#>  [5,]  212  222  213  223  214  224  215  225
#>  [6,]  312  322  313  323  314  324  315  325
#>  [7,]  113  123  114  124  115  125  116  126
#>  [8,]  213  223  214  224  215  225  216  226
#>  [9,]  313  323  314  324  315  325  316  326
bhankel(test_array(dim = c(2,2,6)), d = c(3,7))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]  111  121  112  122  113  123  114  124  115   125   116   126     0     0
#> [2,]  211  221  212  222  213  223  214  224  215   225   216   226     0     0
#> [3,]  112  122  113  123  114  124  115  125  116   126     0     0     0     0
#> [4,]  212  222  213  223  214  224  215  225  216   226     0     0     0     0
#> [5,]  113  123  114  124  115  125  116  126    0     0     0     0     0     0
#> [6,]  213  223  214  224  215  225  216  226    0     0     0     0     0     0
bhankel(test_array(dim = c(1,2,6)), d = c(1,2))
#>      [,1] [,2] [,3] [,4]
#> [1,]  111  121  112  122
bhankel(test_array(dim = c(2,3)), d = c(2,2))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]   11   12   13    0    0    0
#> [2,]   21   22   23    0    0    0
#> [3,]    0    0    0    0    0    0
#> [4,]    0    0    0    0    0    0
bhankel(test_array(dim = c(2,2,6)), d = c(3,0))
#>     
#> [1,]
#> [2,]
#> [3,]
#> [4,]
#> [5,]
#> [6,]
bhankel(test_array(dim = c(2,2,0)), d = c(2,2))
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0    0    0
#> [2,]    0    0    0    0
#> [3,]    0    0    0    0
#> [4,]    0    0    0    0

if (FALSE) { # \dontrun{
# the following examples throw an error
bhankel(1:5, d = c(-1,1))
bhankel(test_array(dim = c(2,3,2,1)))
} # }
```
