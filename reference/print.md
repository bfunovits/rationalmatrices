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
#> [1,]     1.05  1.03  1.29    -0.33  0.17  1.16    -0.45  0.88  0.96
#> [2,]    -2.18 -0.74 -0.67    -0.63 -0.67 -0.78     0.80 -0.14 -0.40
#> 
#> format = i|zj 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>      [,1] z^0   z^1   z^2 [,2] z^0   z^1   z^2 [,3] z^0   z^1   z^2
#> [1,]     1.05 -0.33 -0.45     1.03  0.17  0.88     1.29  1.16  0.96
#> [2,]    -2.18 -0.63  0.80    -0.74 -0.67 -0.14    -0.67 -0.78 -0.40
#> 
#> format = iz|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> z^0 [1,]  1.05  1.03  1.29
#>     [2,] -2.18 -0.74 -0.67
#> z^1 [1,] -0.33  0.17  1.16
#>     [2,] -0.63 -0.67 -0.78
#> z^2 [1,] -0.45  0.88  0.96
#>     [2,]  0.80 -0.14 -0.40
#> 
#> format = zi|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> [1,] z^0  1.05  1.03  1.29
#>      z^1 -0.33  0.17  1.16
#>      z^2 -0.45  0.88  0.96
#> [2,] z^0 -2.18 -0.74 -0.67
#>      z^1 -0.63 -0.67 -0.78
#>      z^2  0.80 -0.14 -0.40
#> 
#> format = i|j|z 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#> , , z^0
#> 
#>       [,1]  [,2]  [,3]
#> [1,]  1.05  1.03  1.29
#> [2,] -2.18 -0.74 -0.67
#> 
#> , , z^1
#> 
#>       [,1]  [,2]  [,3]
#> [1,] -0.33  0.17  1.16
#> [2,] -0.63 -0.67 -0.78
#> 
#> , , z^2
#> 
#>       [,1]  [,2]  [,3]
#> [1,] -0.45  0.88  0.96
#> [2,]  0.80 -0.14 -0.40
#> 
#> 
#> format = character 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>                         [,1]                     [,2]                    [,3]
#> [1,]  1.05 - 0.33z - 0.45z^2   1.03 + 0.17z + 0.88z^2  1.29 + 1.16z + 0.96z^2
#> [2,]  -2.18 - 0.63z + 0.8z^2  -0.74 - 0.67z - 0.14z^2  -0.67 - 0.78z - 0.4z^2 

# "empty" (2 x 0) polynomial matrix (degree = 2)
a = test_polm(dim = c(2,0), degree = 0)
print(a)
#> ( 2 x 0 ) matrix polynomial with degree <= -1 

# random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
a = polm(array(complex(real = stats::rnorm(2*1*3), 
                       imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
print(a, digits = 2)
#> ( 2 x 1 ) matrix polynomial with degree <= 2 
#>         z^0 [,1]    z^1 [,1]   z^2 [,1]
#> [1,] -0.14+1.08i -1.48-0.39i 0.78-0.85i
#> [2,]  0.49-0.89i -1.20-0.46i 0.78-1.05i
if (FALSE) { # \dontrun{
# the format option 'character' is only implemented for polynomials matrices 
# with real coefficients!
print(a, digits = 2, format = 'character')
} # }

# print a rational matrix in statespace form
a = test_stsp(dim = c(3,3), s = 2)
print(a, digits = 2)
#> statespace realization [3,3] with s = 2 states
#>      s[1]  s[2]  u[1] u[2]  u[3]
#> s[1] 0.53  0.64 -1.50 1.71  1.71
#> s[2] 0.63  0.64 -1.09 0.11 -0.29
#> x[1] 1.26  1.43  1.00 0.00  0.00
#> x[2] 1.67 -0.84  0.00 1.00  0.00
#> x[3] 2.25 -1.35  0.00 0.00  1.00

# print a rational matrix in 'lmfd' form 
a = test_lmfd(dim = c(2,3), degrees = c(2,1))
print(a, digits = 2, format = 'character')
#> ( 2 x 3 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 1)
#> left factor a(z):
#>                      [,1]                 [,2]
#> [1,]  1 - 1.51z - 1.18z^2      0.09z + 2.07z^2
#> [2,]     -0.98z + 0.32z^2  1 + 0.88z + 0.26z^2 
#> right factor b(z):
#>                [,1]          [,2]           [,3]
#> [1,]  -1.03 + 0.53z  0.19 + 0.41z  -0.05 - 0.62z
#> [2,]  -0.02 + 0.33z  0.81 - 1.48z    1.09 + 2.5z 

# print impulse response 
print(pseries(a), format = 'i|zj', digits = 2)
#> ( 2 x 3 ) impulse response with maximum lag = 5 
#>      [,1] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5 [,2] lag=0  lag=1  lag=2
#> [1,]      -1.03  -1.03  -2.67  -3.88  -8.64 -13.34       0.19   0.62  -0.34
#> [2,]      -0.02  -0.66  -0.10  -2.03  -1.14  -5.70       0.81  -2.01   2.09
#>       lag=3  lag=4  lag=5 [,3] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5
#> [1,]   4.18   1.74  10.92      -0.05   -0.8  -3.65  -9.35 -13.41 -27.40
#> [2,]  -1.84   5.27  -3.76       1.09    1.5  -2.36  -1.64  -5.94  -4.53

# print frequency response 
print(zvalues(a), format = 'iz|j', digits = 2)
#> ( 2 x 3 ) frequency response
#>                             [,1]        [,2]        [,3]
#>          z=1+0i [1,]  0.79+0.00i -1.24+0.00i  4.19+0.00i
#>                 [2,]  0.39+0.00i -0.70+0.00i  2.98+0.00i
#>  z=0.309-0.951i [1,] -0.35+0.20i  0.06+1.10i  2.22-1.94i
#>                 [2,] -0.01+0.05i -0.32+1.56i  3.24-1.80i
#> z=-0.809-0.588i [1,] -0.37-0.06i  1.18-0.61i -1.00-0.08i
#>                 [2,] -0.25+0.32i  0.26+1.22i  0.46-1.08i
#> z=-0.809+0.588i [1,] -0.37+0.06i  1.18+0.61i -1.00+0.08i
#>                 [2,] -0.25-0.32i  0.26-1.22i  0.46+1.08i
#>  z=0.309+0.951i [1,] -0.35-0.20i  0.06-1.10i  2.22+1.94i
#>                 [2,] -0.01-0.05i -0.32-1.56i  3.24+1.80i
```
