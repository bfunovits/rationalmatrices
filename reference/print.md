# Print Methods

Printing rational matrix objects.

## Usage

``` r
# S3 method for class 'lpolm'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  ...
)

# S3 method for class 'polm'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  ...
)

# S3 method for class 'lmfd'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  ...
)

# S3 method for class 'rmfd'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character"),
  ...
)

# S3 method for class 'stsp'
print(x, digits = NULL, ...)

# S3 method for class 'pseries'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)

# S3 method for class 'zvalues'
print(
  x,
  digits = NULL,
  format = c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z"),
  ...
)
```

## Arguments

- x:

  rational matrix object, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  object.

- digits:

  (integer) if non `NULL` then correspondingly rounded numbers are
  printed, see [`round`](https://rdrr.io/r/base/Round.html).

- format:

  (character string) selects specific output formats. Note that
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  objects have no format option. The option `'character'` is only
  implemented for polynomials, Laurent polynomials, LMFDs and RMFDs with
  real coefficients, se

- ...:

  Further parameters are ignored.

## Value

`invisible(x)`

## Examples

``` r
# for polynomials six different print formats are implemented ###################
a = test_polm(dim = c(2,3), degree = 2, random = TRUE)

for (fmt in c("i|jz", "i|zj", "iz|j", "zi|j", "i|j|z", "character")) {
   cat('\nformat =', fmt, '\n')
   print(a, digits = 2, format = fmt)
}
#> 
#> format = i|jz 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]     2.25 -2.07 -0.28    -0.37  1.88  1.03     1.19 -0.57 -0.38
#> [2,]    -1.23  0.13  0.46    -1.91  0.43  0.03     0.66 -0.57  2.84
#> 
#> format = i|zj 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>      [,1] z^0   z^1  z^2 [,2] z^0  z^1   z^2 [,3] z^0  z^1   z^2
#> [1,]     2.25 -0.37 1.19    -2.07 1.88 -0.57    -0.28 1.03 -0.38
#> [2,]    -1.23 -1.91 0.66     0.13 0.43 -0.57     0.46 0.03  2.84
#> 
#> format = iz|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> z^0 [1,]  2.25 -2.07 -0.28
#>     [2,] -1.23  0.13  0.46
#> z^1 [1,] -0.37  1.88  1.03
#>     [2,] -1.91  0.43  0.03
#> z^2 [1,]  1.19 -0.57 -0.38
#>     [2,]  0.66 -0.57  2.84
#> 
#> format = zi|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> [1,] z^0  2.25 -2.07 -0.28
#>      z^1 -0.37  1.88  1.03
#>      z^2  1.19 -0.57 -0.38
#> [2,] z^0 -1.23  0.13  0.46
#>      z^1 -1.91  0.43  0.03
#>      z^2  0.66 -0.57  2.84
#> 
#> format = i|j|z 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#> , , z^0
#> 
#>       [,1]  [,2]  [,3]
#> [1,]  2.25 -2.07 -0.28
#> [2,] -1.23  0.13  0.46
#> 
#> , , z^1
#> 
#>       [,1] [,2] [,3]
#> [1,] -0.37 1.88 1.03
#> [2,] -1.91 0.43 0.03
#> 
#> , , z^2
#> 
#>      [,1]  [,2]  [,3]
#> [1,] 1.19 -0.57 -0.38
#> [2,] 0.66 -0.57  2.84
#> 
#> 
#> format = character 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>                          [,1]                     [,2]                     [,3]
#> [1,]   2.25 - 0.37z + 1.19z^2  -2.07 + 1.88z - 0.57z^2  -0.28 + 1.03z - 0.38z^2
#> [2,]  -1.23 - 1.91z + 0.66z^2   0.13 + 0.43z - 0.57z^2   0.46 + 0.03z + 2.84z^2 

# "empty" (2 x 0) polynomial matrix (degree = 2)
a = test_polm(dim = c(2,0), degree = 0)
print(a)
#> ( 2 x 0 ) matrix polynomial with degree <= -1 

# random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
a = polm(array(complex(real = stats::rnorm(2*1*3), 
                       imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
print(a, digits = 2)
#> ( 2 x 1 ) matrix polynomial with degree <= 2 
#>        z^0 [,1]    z^1 [,1]   z^2 [,1]
#> [1,] 0.68+0.62i -0.57-0.26i 0.93+0.30i
#> [2,] 1.26+0.05i  0.26-0.67i 0.87+1.25i
if (FALSE) { # \dontrun{
# the format option 'character' is only implemented for polynomials matrices 
# with real coefficients!
print(a, digits = 2, format = 'character')
} # }

# print a rational matrix in statespace form
a = test_stsp(dim = c(3,3), s = 2)
print(a, digits = 2)
#> statespace realization [3,3] with s = 2 states
#>       s[1]  s[2]  u[1]  u[2]  u[3]
#> s[1] -0.07 -0.73 -1.48  1.51 -0.29
#> s[2] -0.08  0.75  1.19 -1.44 -1.39
#> x[1]  0.10  0.12  1.00  0.00  0.00
#> x[2]  0.48  0.45  0.00  1.00  0.00
#> x[3] -0.56 -0.08  0.00  0.00  1.00

# print a rational matrix in 'lmfd' form 
a = test_lmfd(dim = c(2,3), degrees = c(2,1))
print(a, digits = 2, format = 'character')
#> ( 2 x 3 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 1)
#> left factor a(z):
#>                      [,1]                 [,2]
#> [1,]  1 - 1.08z + 1.12z^2      0.52z + 1.04z^2
#> [2,]      0.36z + 0.13z^2  1 - 0.27z + 0.99z^2 
#> right factor b(z):
#>                [,1]           [,2]           [,3]
#> [1,]  -1.38 - 0.32z   0.47 + 0.37z  -0.22 + 0.06z
#> [2,]   1.21 + 1.09z  -1.25 - 0.68z  -0.44 + 0.62z 

# print impulse response 
print(pseries(a), format = 'i|zj', digits = 2)
#> ( 2 x 3 ) impulse response with maximum lag = 5 
#>      [,1] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5 [,2] lag=0  lag=1  lag=2
#> [1,]      -1.38  -2.43  -3.31  -3.03   0.23   3.40       0.47   1.52   3.02
#> [2,]       1.21   1.91   0.36  -0.31   1.05   0.89      -1.25  -1.18   0.31
#>       lag=3  lag=4  lag=5 [,3] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5
#> [1,]   2.62  -0.87  -3.02      -0.22   0.06   0.47  -0.47  -1.36   0.01
#> [2,]  -0.01  -1.62  -0.45      -0.44   0.58   0.60  -0.58  -0.64   0.94

# print frequency response 
print(zvalues(a), format = 'iz|j', digits = 2)
#> ( 2 x 3 ) frequency response
#>                             [,1]        [,2]        [,3]
#>          z=1+0i [1,] -6.25+0.00i  4.27+0.00i -0.52+0.00i
#>                 [2,]  3.09+0.00i -2.32+0.00i  0.25+0.00i
#>  z=0.309-0.951i [1,]  2.87+3.53i -2.29-3.17i  1.02-0.84i
#>                 [2,] -0.26+0.29i  0.63-0.63i  0.52+0.08i
#> z=-0.809-0.588i [1,] -0.34+0.36i  0.10-0.11i  0.01+0.14i
#>                 [2,] -0.09-0.32i -0.17+0.37i -0.51+0.15i
#> z=-0.809+0.588i [1,] -0.34-0.36i  0.10+0.11i  0.01-0.14i
#>                 [2,] -0.09+0.32i -0.17-0.37i -0.51-0.15i
#>  z=0.309+0.951i [1,]  2.87-3.53i -2.29+3.17i  1.02+0.84i
#>                 [2,] -0.26-0.29i  0.63+0.63i  0.52-0.08i
```
