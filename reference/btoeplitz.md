# Block Toeplitz matrix

Construct a block Toeplitz matrix from two 3-dimensional arrays `R` and
`C`. The array `R` determines the first block row and `C` the first
block column of the result. The \\(i,j)\\-th block of the Toeplitz
matrix is `R[,,j-i+1]` for \\j\geq i\\ and `C[,,i-j+1]` for \\j \< i\\.
In particular, note that the \\(1,1)\\ block is set to `R[,,1]` (while
`C[,,1]` is ignored).

## Usage

``` r
btoeplitz(R, C)
```

## Arguments

- R:

  Array with dimensions \\(p,q,n)\\. Corresponds to the first (block-)
  row, containing \\n\\ matrices of dimension \\(p \times q)\\. A vector
  is coerced to an \\(1,1,n)\\ dimensional array and a \\(p,q)\\ matrix
  is interpreted as an \\(p,q,1)\\ dimensional array.

- C:

  Array with dimensions \\(p,q,m)\\. Corresponds to the first (block-)
  column, containing \\m\\ matrices of dimension \\(p \times q)\\. A
  vector is coerced to an \\(1,1,m)\\ dimensional array and a \\(p,q)\\
  matrix is interpreted as an \\(p,q,1)\\ dimensional array.

## Value

Matrix of size \\( pm \times qn )\\.

## Details

If only the argument `R` is provided then `C` is set to
`C = aperm(R,c(2,1,3)))`. So `btoeplitz(R = X)` is equivalent to
`btoeplitz(R = X, C = aperm(X,c(2,1,3)))`. Analogously
`btoeplitz(C = X)` is equivalent to
`btoeplitz(R = aperm(X,c(2,1,3)), C = X)`. In these cases \\p=q\\ must
hold. The so constructed Toeplitz matrix is symmetric if and only if
`X[,,1]` is symmetric. Note also that even in the case where only `C` is
supplied the "`R` wins" rule holds, i.e. `btoeplitz(C = X)` sets the
\\(1,1)\\ block equal to `t(X[,,1])`.

## Examples

``` r
btoeplitz(0:3)
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    2    3
#> [2,]    1    0    1    2
#> [3,]    2    1    0    1
#> [4,]    3    2    1    0
btoeplitz(0:3,-(0:3))
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    1    2    3
#> [2,]   -1    0    1    2
#> [3,]   -2   -1    0    1
#> [4,]   -3   -2   -1    0

btoeplitz(R = test_array(dim=c(2,3,3)), C = -test_array(dim = c(2,3,2)))
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
#> [1,]  111  121  131  112  122  132  113  123  133
#> [2,]  211  221  231  212  222  232  213  223  233
#> [3,] -112 -122 -132  111  121  131  112  122  132
#> [4,] -212 -222 -232  211  221  231  212  222  232
btoeplitz(R = test_array(dim=c(2,2,1)), C = test_array(dim=c(2,2,0)))
#>      [,1] [,2]
btoeplitz(R = test_array(dim=c(2,2,3)))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  111  121  112  122  113  123
#> [2,]  211  221  212  222  213  223
#> [3,]  112  212  111  121  112  122
#> [4,]  122  222  211  221  212  222
#> [5,]  113  213  112  212  111  121
#> [6,]  123  223  122  222  211  221
btoeplitz(C = -test_array(dim=c(2,2,3)))
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,] -111 -211 -112 -212 -113 -213
#> [2,] -121 -221 -122 -222 -123 -223
#> [3,] -112 -122 -111 -211 -112 -212
#> [4,] -212 -222 -121 -221 -122 -222
#> [5,] -113 -123 -112 -122 -111 -211
#> [6,] -213 -223 -212 -222 -121 -221
# create a symmetric matrix
X = test_array(dim= c(2,2,3))
X[,,1] = (X[,,1] + t(X[,,1]))/2
btoeplitz(R = X)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  111  166  112  122  113  123
#> [2,]  166  221  212  222  213  223
#> [3,]  112  212  111  166  112  122
#> [4,]  122  222  166  221  212  222
#> [5,]  113  213  112  212  111  166
#> [6,]  123  223  122  222  166  221
btoeplitz(C = X)
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]  111  166  112  212  113  213
#> [2,]  166  221  122  222  123  223
#> [3,]  112  122  111  166  112  212
#> [4,]  212  222  166  221  122  222
#> [5,]  113  123  112  122  111  166
#> [6,]  213  223  212  222  166  221

if (FALSE) { # \dontrun{
# the following examples throw an error
btoeplitz(R = test_array(dim=c(2,1,3)), C = -test_array(dim = c(2,2,2)))
btoeplitz(R = test_array(dim=c(2,1,3)))
btoeplitz(C = test_array(dim=c(2,1,3)))
} # }
```
