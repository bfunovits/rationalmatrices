# Replace Parts of a (Laurent) Polynomial Matrix

The assignment operation `x[,] <- value` for (Laurent) polynomial
matrices works quite analogously to the assignment operation of
"ordinary" matrices.  
In the case of Laurent polynomial objects
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
and if `value` is not an
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object itself (i.e. `value` is a vector, matrix, or array), `value` will
first be coerced to an
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object with `min_deg = 0`.

## Usage

``` r
# S3 method for class 'polm'
x[i, j] <- value

# S3 method for class 'lpolm'
x[i, j] <- value
```

## Arguments

- x:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  or
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  object

- i, j:

  indices

- value:

  Either a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  or
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  object, or a vector/matrix/array which may be coerced to a `polm` or
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  object by `polm(value)` or `lpolm(value)`.

## Details

Note that "named" arguments are not supported (in order to simplify the
coding).

## Examples

``` r
a = test_polm(dim = c(3,2), degree = 1)
print(a, format = 'c')
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>             [,1]        [,2]
#> [1,]  110 + 111z  120 + 121z
#> [2,]  210 + 211z  220 + 221z
#> [3,]  310 + 311z  320 + 321z 

a[FALSE] = 0           # no items to replace,
print(a, format = 'c') # a is not changed
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>             [,1]        [,2]
#> [1,]  110 + 111z  120 + 121z
#> [2,]  210 + 211z  220 + 221z
#> [3,]  310 + 311z  320 + 321z 

a[lower.tri(matrix(0, nrow = 3, ncol = 2))] = 0 # set elements below the 
print(a, format = 'c')                          # diagonal equal to zero
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>             [,1]        [,2]
#> [1,]  110 + 111z  120 + 121z
#> [2,]           0  220 + 221z
#> [3,]           0           0 

a[3,1] = c(1,-1)       # set (3,1) element
print(a, format = 'c') # equal to (1 - z)
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>             [,1]        [,2]
#> [1,]  110 + 111z  120 + 121z
#> [2,]           0  220 + 221z
#> [3,]       1 - z           0 

a[1:2, 2:1] = c(0,1)   # set the elements in rows 1,2 and coluimns 1,2
print(a, format = 'c') # equal to z
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>        [,1]  [,2]
#> [1,]      z     z
#> [2,]      z     z
#> [3,]  1 - z     0 

a[2, ] = test_polm(dim = c(1,2), degree = 4)
print(a, format = 'c')
#> ( 3 x 2 ) matrix polynomial with degree <= 4 
#>                                        [,1]                                   [,2]
#> [1,]                                      z                                      z
#> [2,]  110 + 111z + 112z^2 + 113z^3 + 114z^4  120 + 121z + 122z^2 + 123z^3 + 124z^4
#> [3,]                                  1 - z                                      0 

a[, 1] = test_polm(dim = c(2,1), degree = 4) # this gives a warning
#> Warning: number of items to replace is not a multiple of replacement length
print(a, format = 'c')
#> ( 3 x 2 ) matrix polynomial with degree <= 4 
#>                                        [,1]                                   [,2]
#> [1,]  110 + 111z + 112z^2 + 113z^3 + 114z^4                                      z
#> [2,]  210 + 211z + 212z^2 + 213z^3 + 214z^4  120 + 121z + 122z^2 + 123z^3 + 124z^4
#> [3,]  110 + 111z + 112z^2 + 113z^3 + 114z^4                                      0 

if (FALSE) { # \dontrun{
a[i=1] = 0   # named arguments are not supported
} # }  
(lp = lpolm(1:3, min_deg = -1))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1]
#> [1,]         1        2        3
lp[1] = 0
lp
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1]
#> [1,]         0        0        0

(lp = lpolm(1:3, min_deg = -1))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1]
#> [1,]         1        2        3
lp[1] = polm(1:4)
lp
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 3, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1]
#> [1,]         0        1        2        3        4

(lp = lpolm(1:3, min_deg = -1))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1]
#> [1,]         1        2        3
lp[1] = lpolm(1, min_deg = 0)
lp
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1]
#> [1,]         0        1        0
```
