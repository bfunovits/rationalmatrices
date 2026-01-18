# Companion Matrix of a Polynomial Matrix

Computes a companion matrix for a square (\\m,m)\\-dimensional), matrix
polynomial \$\$a(z) = a_0 + a_1 z + \cdots + a_p z^p\$\$ The companion
matrix is e.g. used to determine the zeroes of a polynomial matrix, see
[`zeroes`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md).  
Note that the function throws an error, if the constant term \\a_0\\ is
singular. There is no check whether some of the leading coefficients are
zero. So the results is an \\(mp,mp)\\-dimensional matrix, even if
\\a_p\\ is zero.

## Usage

``` r
companion_matrix(a)
```

## Arguments

- a:

  A square polynomial matrix, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

## Value

A (companion) matrix of dimensions \\(mp,mp)\\.

## Examples

``` r
companion_matrix(polm(c(1,0,0,0.5,0))) # scalar polynomial
#>      [,1] [,2] [,3] [,4]
#> [1,]    0    0 -0.5    0
#> [2,]    1    0  0.0    0
#> [3,]    0    1  0.0    0
#> [4,]    0    0  1.0    0
companion_matrix(polm(diag(3)))        # zero degree polynomial 
#> <0 x 0 matrix>
companion_matrix(polm(dbind(d = 3, diag(2), -test_array(dim = c(2,2,1)))))
#>      [,1] [,2]
#> [1,]  111  121
#> [2,]  211  221
companion_matrix(polm(dbind(d = 3, diag(2), -test_array(dim = c(2,2,2)))))
#>      [,1] [,2] [,3] [,4]
#> [1,]  111  121  112  122
#> [2,]  211  221  212  222
#> [3,]    1    0    0    0
#> [4,]    0    1    0    0

if (FALSE) { # \dontrun{
# the following examples throw an error
companion_matrix(polm(c(0,0,0,0.5))) # constant term is zero
companion_matrix(polm(test_array(dim = c(2,1,3)))) # non-square polynomial
companion_matrix(polm(test_array(dim = c(2,2,0)))) # zero polynomial
} # }
```
