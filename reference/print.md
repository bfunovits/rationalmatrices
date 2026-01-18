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
#> [1,]     0.15 -0.13  0.57     0.06 -1.86  1.29    -1.06  0.57  0.56
#> [2,]    -1.49 -0.37 -0.35     1.72  0.65  0.66     2.28  0.33  1.35
#> 
#> format = i|zj 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>      [,1] z^0  z^1   z^2 [,2] z^0   z^1  z^2 [,3] z^0  z^1  z^2
#> [1,]     0.15 0.06 -1.06    -0.13 -1.86 0.57     0.57 1.29 0.56
#> [2,]    -1.49 1.72  2.28    -0.37  0.65 0.33    -0.35 0.66 1.35
#> 
#> format = iz|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> z^0 [1,]  0.15 -0.13  0.57
#>     [2,] -1.49 -0.37 -0.35
#> z^1 [1,]  0.06 -1.86  1.29
#>     [2,]  1.72  0.65  0.66
#> z^2 [1,] -1.06  0.57  0.56
#>     [2,]  2.28  0.33  1.35
#> 
#> format = zi|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> [1,] z^0  0.15 -0.13  0.57
#>      z^1  0.06 -1.86  1.29
#>      z^2 -1.06  0.57  0.56
#> [2,] z^0 -1.49 -0.37 -0.35
#>      z^1  1.72  0.65  0.66
#>      z^2  2.28  0.33  1.35
#> 
#> format = i|j|z 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#> , , z^0
#> 
#>       [,1]  [,2]  [,3]
#> [1,]  0.15 -0.13  0.57
#> [2,] -1.49 -0.37 -0.35
#> 
#> , , z^1
#> 
#>      [,1]  [,2] [,3]
#> [1,] 0.06 -1.86 1.29
#> [2,] 1.72  0.65 0.66
#> 
#> , , z^2
#> 
#>       [,1] [,2] [,3]
#> [1,] -1.06 0.57 0.56
#> [2,]  2.28 0.33 1.35
#> 
#> 
#> format = character 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>                          [,1]                     [,2]                     [,3]
#> [1,]   0.15 + 0.06z - 1.06z^2  -0.13 - 1.86z + 0.57z^2   0.57 + 1.29z + 0.56z^2
#> [2,]  -1.49 + 1.72z + 2.28z^2  -0.37 + 0.65z + 0.33z^2  -0.35 + 0.66z + 1.35z^2 

# "empty" (2 x 0) polynomial matrix (degree = 2)
a = test_polm(dim = c(2,0), degree = 0)
print(a)
#> ( 2 x 0 ) matrix polynomial with degree <= -1 

# random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
a = polm(array(complex(real = stats::rnorm(2*1*3), 
                       imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
print(a, digits = 2)
#> ( 2 x 1 ) matrix polynomial with degree <= 2 
#>         z^0 [,1]   z^1 [,1]    z^2 [,1]
#> [1,] -0.34+0.27i 0.41+0.00i  0.56-1.03i
#> [2,] -0.43-0.90i 2.02+1.71i -0.10+0.14i
if (FALSE) { # \dontrun{
# the format option 'character' is only implemented for polynomials matrices 
# with real coefficients!
print(a, digits = 2, format = 'character')
} # }

# print a rational matrix in statespace form
a = test_stsp(dim = c(3,3), s = 2)
print(a, digits = 2)
#> statespace realization [3,3] with s = 2 states
#>       s[1]  s[2] u[1] u[2]  u[3]
#> s[1] -0.19 -0.52 1.17 1.32 -1.12
#> s[2] -0.09  1.46 0.30 0.03 -0.06
#> x[1]  0.49 -0.64 1.00 0.00  0.00
#> x[2]  0.28 -0.15 0.00 1.00  0.00
#> x[3]  1.89  0.39 0.00 0.00  1.00

# print a rational matrix in 'lmfd' form 
a = test_lmfd(dim = c(2,3), degrees = c(2,1))
print(a, digits = 2, format = 'character')
#> ( 2 x 3 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 1)
#> left factor a(z):
#>                      [,1]                 [,2]
#> [1,]  1 + 1.16z + 0.35z^2      0.59z - 0.66z^2
#> [2,]     -0.48z + 1.24z^2  1 + 0.18z + 0.73z^2 
#> right factor b(z):
#>               [,1]           [,2]          [,3]
#> [1,]   1.66 + 0.7z   0.24 - 0.58z  1.06 + 0.47z
#> [2,]  1.72 - 0.71z  -1.32 - 0.21z  -1.16 - 0.3z 

# print impulse response 
print(pseries(a), format = 'i|zj', digits = 2)
#> ( 2 x 3 ) impulse response with maximum lag = 5 
#>      [,1] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5 [,2] lag=0  lag=1  lag=2
#> [1,]       1.66  -2.24   3.28  -0.58  -6.49  12.44       0.24  -0.08  -0.95
#> [2,]       1.72  -0.21  -4.36   5.29  -2.07  -5.95      -1.32   0.14   0.61
#>       lag=3  lag=4  lag=5 [,3] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5
#> [1,]   0.85   0.09  -1.51       1.06  -0.07  -1.30   2.13  -1.95  -0.87
#> [2,]  -0.57   1.24  -0.81      -1.16   0.41  -0.56  -0.75   3.19  -3.59

# print frequency response 
print(zvalues(a), format = 'iz|j', digits = 2)
#> ( 2 x 3 ) frequency response
#>                             [,1]        [,2]        [,3]
#>          z=1+0i [1,]  0.94+0.00i -0.16+0.00i  0.58+0.00i
#>                 [2,]  0.16+0.00i -0.74+0.00i -0.99+0.00i
#>  z=0.309-0.951i [1,]  0.86-0.98i  0.51+0.82i  1.23+0.33i
#>                 [2,]  2.45+2.71i -2.16-0.04i -0.98+0.88i
#> z=-0.809-0.588i [1,]  1.12-2.14i -0.08+0.26i -0.43+0.08i
#>                 [2,] -1.03+0.97i -0.49+0.31i -0.04+0.71i
#> z=-0.809+0.588i [1,]  1.12+2.14i -0.08-0.26i -0.43-0.08i
#>                 [2,] -1.03-0.97i -0.49-0.31i -0.04-0.71i
#>  z=0.309+0.951i [1,]  0.86+0.98i  0.51-0.82i  1.23-0.33i
#>                 [2,]  2.45-2.71i -2.16+0.04i -0.98-0.88i
```
