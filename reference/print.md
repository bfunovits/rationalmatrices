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
#> [1,]    -0.83 -1.13  0.92    -0.72 -1.28  0.43     0.65  0.46 -0.17
#> [2,]     0.61  0.72  0.51     0.44  0.69  1.72    -1.35  0.83  1.50
#> 
#> format = i|zj 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>      [,1] z^0   z^1   z^2 [,2] z^0   z^1  z^2 [,3] z^0  z^1   z^2
#> [1,]    -0.83 -0.72  0.65    -1.13 -1.28 0.46     0.92 0.43 -0.17
#> [2,]     0.61  0.44 -1.35     0.72  0.69 0.83     0.51 1.72  1.50
#> 
#> format = iz|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> z^0 [1,] -0.83 -1.13  0.92
#>     [2,]  0.61  0.72  0.51
#> z^1 [1,] -0.72 -1.28  0.43
#>     [2,]  0.44  0.69  1.72
#> z^2 [1,]  0.65  0.46 -0.17
#>     [2,] -1.35  0.83  1.50
#> 
#> format = zi|j 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> [1,] z^0 -0.83 -1.13  0.92
#>      z^1 -0.72 -1.28  0.43
#>      z^2  0.65  0.46 -0.17
#> [2,] z^0  0.61  0.72  0.51
#>      z^1  0.44  0.69  1.72
#>      z^2 -1.35  0.83  1.50
#> 
#> format = i|j|z 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#> , , z^0
#> 
#>       [,1]  [,2] [,3]
#> [1,] -0.83 -1.13 0.92
#> [2,]  0.61  0.72 0.51
#> 
#> , , z^1
#> 
#>       [,1]  [,2] [,3]
#> [1,] -0.72 -1.28 0.43
#> [2,]  0.44  0.69 1.72
#> 
#> , , z^2
#> 
#>       [,1] [,2]  [,3]
#> [1,]  0.65 0.46 -0.17
#> [2,] -1.35 0.83  1.50
#> 
#> 
#> format = character 
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>                          [,1]                     [,2]                    [,3]
#> [1,]  -0.83 - 0.72z + 0.65z^2  -1.13 - 1.28z + 0.46z^2  0.92 + 0.43z - 0.17z^2
#> [2,]   0.61 + 0.44z - 1.35z^2   0.72 + 0.69z + 0.83z^2   0.51 + 1.72z + 1.5z^2 

# "empty" (2 x 0) polynomial matrix (degree = 2)
a = test_polm(dim = c(2,0), degree = 0)
print(a)
#> ( 2 x 0 ) matrix polynomial with degree <= -1 

# random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
a = polm(array(complex(real = stats::rnorm(2*1*3), 
                       imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
print(a, digits = 2)
#> ( 2 x 1 ) matrix polynomial with degree <= 2 
#>         z^0 [,1]    z^1 [,1]    z^2 [,1]
#> [1,] -0.34-0.70i -0.69-1.45i  0.09-1.47i
#> [2,]  1.25+1.53i  0.82-0.74i -0.30-0.75i
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
#> s[1] -0.55  0.38 -0.05 -0.34 -2.20
#> s[2] -0.51 -0.14  2.46 -0.70 -0.23
#> x[1]  0.04 -0.11  1.00  0.00  0.00
#> x[2]  1.45  0.54  0.00  1.00  0.00
#> x[3] -0.12  0.37  0.00  0.00  1.00

# print a rational matrix in 'lmfd' form 
a = test_lmfd(dim = c(2,3), degrees = c(2,1))
print(a, digits = 2, format = 'character')
#> ( 2 x 3 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 1)
#> left factor a(z):
#>                      [,1]                 [,2]
#> [1,]  1 + 0.26z - 1.09z^2     -1.63z + 2.33z^2
#> [2,]     -0.56z + 0.65z^2  1 - 0.14z + 0.53z^2 
#> right factor b(z):
#>                [,1]          [,2]           [,3]
#> [1,]   0.63 + 0.42z  0.72 + 0.83z   -0.82 - 1.4z
#> [2,]  -0.72 - 1.12z  1.23 + 1.35z  -0.91 + 0.48z 

# print impulse response 
print(pseries(a), format = 'i|zj', digits = 2)
#> ( 2 x 3 ) impulse response with maximum lag = 5 
#>      [,1] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5 [,2] lag=0  lag=1  lag=2
#> [1,]       0.63  -0.92   1.19  -0.37   5.61  -6.35       0.72   2.65   0.38
#> [2,]      -0.72  -0.86  -0.66   1.64  -0.41   2.45       1.23   1.92   0.62
#>       lag=3  lag=4  lag=5 [,3] lag=0  lag=1  lag=2  lag=3  lag=4  lag=5
#> [1,]  -0.68  -4.86   4.16      -0.82  -2.67   1.74  -3.92   8.45 -17.19
#> [2,]  -2.46  -1.29  -1.14      -0.91  -0.10  -0.48   2.70  -2.70   5.48

# print frequency response 
print(zvalues(a), format = 'iz|j', digits = 2)
#> ( 2 x 3 ) frequency response
#>                             [,1]        [,2]         [,3]
#>          z=1+0i [1,] 16.76+0.00i  2.21+0.00i -17.12+0.00i
#>                 [2,] -2.44+0.00i  1.69+0.00i   0.85+0.00i
#>  z=0.309-0.951i [1,]  4.68-1.59i -8.72+2.04i   3.10+1.91i
#>                 [2,]  3.79-0.08i -7.91-0.04i   2.63+1.73i
#> z=-0.809-0.588i [1,]  0.48+0.15i -0.28-0.30i  -0.72+0.72i
#>                 [2,]  0.04+0.06i -0.08-0.22i  -0.11+0.00i
#> z=-0.809+0.588i [1,]  0.48-0.15i -0.28+0.30i  -0.72-0.72i
#>                 [2,]  0.04-0.06i -0.08+0.22i  -0.11+0.00i
#>  z=0.309+0.951i [1,]  4.68+1.59i -8.72-2.04i   3.10-1.91i
#>                 [2,]  3.79+0.08i -7.91+0.04i   2.63-1.73i
```
