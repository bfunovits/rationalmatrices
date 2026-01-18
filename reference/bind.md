# Combine Rational Matrices by Rows or Columns

The function \`rbind()\` and \`cbind()\` take a sequence of rational
matrix objects (i.e.
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
or
[`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
objects) and combines these matrices by rows or columns.

## Usage

``` r
# S3 method for class 'ratm'
rbind(...)

# S3 method for class 'ratm'
cbind(...)
```

## Arguments

- ...:

  rational matrix objects, or objects which may be coerced to rational
  matrix objects.

## Value

rational matrix object.

## Details

The input matrices are first coerced to objects of the same class, as
described for the group operator methods in
[`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md).

## Examples

``` r
a1 = test_polm(dim = c(2,3), degree = 2)
a2 = test_polm(dim = c(1,3), degree = 1)

rbind(a1, a2)                  # => polm object
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]      110   120   130      111   121   131      112   122   132
#> [2,]      210   220   230      211   221   231      212   222   232
#> [3,]      110   120   130      111   121   131        0     0     0
rbind(lmfd(diag(2), a1), a2)   # => stsp object
#> statespace realization [3,3] with s = 5 states
#>      s[1] s[2] s[3] s[4] s[5] u[1] u[2] u[3]
#> s[1]    0    0    0    0    0  112  122  132
#> s[2]    0    0    0    0    0  212  222  232
#> s[3]    1    0    0    0    0  111  121  131
#> s[4]    0    1    0    0    0  211  221  231
#> s[5]    0    0    0    0    0  111  121  131
#> x[1]    0    0    1    0    0  110  120  130
#> x[2]    0    0    0    1    0  210  220  230
#> x[3]    0    0    0    0    1  110  120  130
rbind(a1, as.stsp(a2))         # => stsp object
#> statespace realization [3,3] with s = 5 states
#>      s[1] s[2] s[3] s[4] s[5] u[1] u[2] u[3]
#> s[1]    0    0    1    0    0  111  121  131
#> s[2]    0    0    0    1    0  211  221  231
#> s[3]    0    0    0    0    0  112  122  132
#> s[4]    0    0    0    0    0  212  222  232
#> s[5]    0    0    0    0    0  111  121  131
#> x[1]    1    0    0    0    0  110  120  130
#> x[2]    0    1    0    0    0  210  220  230
#> x[3]    0    0    0    0    1  110  120  130
rbind(a1, pseries(a2))         # pseries object 
#> ( 3 x 3 ) impulse response with maximum lag = 5 
#>      lag=0 [,1]  [,2]  [,3] lag=1 [,1]  [,2]  [,3] lag=2 [,1]  [,2]  [,3]
#> [1,]        110   120   130        111   121   131        112   122   132
#> [2,]        210   220   230        211   221   231        212   222   232
#> [3,]        110   120   130        111   121   131          0     0     0
#>      lag=3 [,1]  [,2]  [,3] lag=4 [,1]  [,2]  [,3] lag=5 [,1]  [,2]  [,3]
#> [1,]          0     0     0          0     0     0          0     0     0
#> [2,]          0     0     0          0     0     0          0     0     0
#> [3,]          0     0     0          0     0     0          0     0     0
rbind(zvalues(a1), a2)        # zvalues object
#> ( 3 x 3 ) frequency response
#>      z[1] [,1]   [,2]   [,3]           z[2] [,1]                [,2]
#> [1,]    333+0i 363+0i 393+0i  53.69098-171.3992i  58.69098-186.7876i
#> [2,]    633+0i 663+0i 693+0i 103.69098-325.2834i 108.69098-340.6718i
#> [3,]    221+0i 241+0i 261+0i 144.30089-105.5673i 157.39106-115.0778i
#>                     [,3]           z[3] [,1]                [,2]
#> [1,]  63.69098-202.1761i  54.80902+41.27417i  59.80902+44.90688i
#> [2,] 113.69098-356.0602i 104.80902+77.60129i 109.80902+81.23401i
#> [3,] 170.48123-124.5884i  20.19911-65.24416i  22.10894-71.12202i
#>                     [,3]           z[4] [,1]                [,2]
#> [1,]  64.80902+48.53959i  54.80902-41.27417i  59.80902-44.90688i
#> [2,] 114.80902+84.86672i 104.80902-77.60129i 109.80902-81.23401i
#> [3,]  24.01877-76.99987i  20.19911+65.24416i  22.10894+71.12202i
#>                     [,3]           z[5] [,1]                [,2]
#> [1,]  64.80902-48.53959i  53.69098+171.3992i  58.69098+186.7876i
#> [2,] 114.80902-84.86672i 103.69098+325.2834i 108.69098+340.6718i
#> [3,]  24.01877+76.99987i 144.30089+105.5673i 157.39106+115.0778i
#>                     [,3]
#> [1,]  63.69098+202.1761i
#> [2,] 113.69098+356.0602i
#> [3,] 170.48123+124.5884i

a2 = test_polm(dim = c(2,1), degree = 3)

cbind(a1, a2)
#> ( 2 x 4 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2]  [,3]  [,4] z^1 [,1]  [,2]  [,3]  [,4] z^2 [,1]  [,2]  [,3]
#> [1,]      110   120   130   110      111   121   131   111      112   122   132
#> [2,]      210   220   230   210      211   221   231   211      212   222   232
#>       [,4] z^3 [,1]  [,2]  [,3]  [,4]
#> [1,]   112        0     0     0   113
#> [2,]   212        0     0     0   213
cbind(lmfd(diag(2), a1), a2)
#> statespace realization [2,4] with s = 7 states
#>      s[1] s[2] s[3] s[4] s[5] s[6] s[7] u[1] u[2] u[3] u[4]
#> s[1]    0    0    0    0    0    0    0  112  122  132    0
#> s[2]    0    0    0    0    0    0    0  212  222  232    0
#> s[3]    1    0    0    0    0    0    0  111  121  131    0
#> s[4]    0    1    0    0    0    0    0  211  221  231    0
#> s[5]    0    0    0    0    0    0    0    0    0    0    1
#> s[6]    0    0    0    0    1    0    0    0    0    0    0
#> s[7]    0    0    0    0    0    1    0    0    0    0    0
#> x[1]    0    0    1    0  111  112  113  110  120  130  110
#> x[2]    0    0    0    1  211  212  213  210  220  230  210
cbind(a1, as.stsp(a2))
#> statespace realization [2,4] with s = 7 states
#>      s[1] s[2] s[3] s[4] s[5] s[6] s[7] u[1] u[2] u[3] u[4]
#> s[1]    0    0    1    0    0    0    0  111  121  131    0
#> s[2]    0    0    0    1    0    0    0  211  221  231    0
#> s[3]    0    0    0    0    0    0    0  112  122  132    0
#> s[4]    0    0    0    0    0    0    0  212  222  232    0
#> s[5]    0    0    0    0    0    0    0    0    0    0    1
#> s[6]    0    0    0    0    1    0    0    0    0    0    0
#> s[7]    0    0    0    0    0    1    0    0    0    0    0
#> x[1]    1    0    0    0  111  112  113  110  120  130  110
#> x[2]    0    1    0    0  211  212  213  210  220  230  210
cbind(a1, pseries(a2))
#> ( 2 x 4 ) impulse response with maximum lag = 5 
#>      lag=0 [,1]  [,2]  [,3]  [,4] lag=1 [,1]  [,2]  [,3]  [,4] lag=2 [,1]  [,2]
#> [1,]        110   120   130   110        111   121   131   111        112   122
#> [2,]        210   220   230   210        211   221   231   211        212   222
#>       [,3]  [,4] lag=3 [,1]  [,2]  [,3]  [,4] lag=4 [,1]  [,2]  [,3]  [,4]
#> [1,]   132   112          0     0     0   113          0     0     0     0
#> [2,]   232   212          0     0     0   213          0     0     0     0
#>      lag=5 [,1]  [,2]  [,3]  [,4]
#> [1,]          0     0     0     0
#> [2,]          0     0     0     0
cbind(zvalues(a1), a2)
#> ( 2 x 4 ) frequency response
#>      z[1] [,1]   [,2]   [,3]   [,4]           z[2] [,1]                [,2]
#> [1,]    333+0i 363+0i 393+0i 446+0i  53.69098-171.3992i  58.69098-186.7876i
#> [2,]    633+0i 663+0i 693+0i 846+0i 103.69098-325.2834i 108.69098-340.6718i
#>                     [,3]                [,4]           z[3] [,1]
#> [1,]  63.69098-202.1761i -37.72794-104.9795i  54.80902+41.27417i
#> [2,] 113.69098-356.0602i -68.62964-200.0851i 104.80902+77.60129i
#>                     [,2]                [,3]                 [,4]
#> [1,]  59.80902+44.90688i  64.80902+48.53959i  89.72794- 66.19522i
#> [2,] 109.80902+81.23401i 114.80902+84.86672i 170.62964-124.97374i
#>                z[4] [,1]                [,2]                [,3]
#> [1,]  54.80902-41.27417i  59.80902-44.90688i  64.80902-48.53959i
#> [2,] 104.80902-77.60129i 109.80902-81.23401i 114.80902-84.86672i
#>                      [,4]           z[5] [,1]                [,2]
#> [1,]  89.72794+ 66.19522i  53.69098+171.3992i  58.69098+186.7876i
#> [2,] 170.62964+124.97374i 103.69098+325.2834i 108.69098+340.6718i
#>                     [,3]                [,4]
#> [1,]  63.69098+202.1761i -37.72794+104.9795i
#> [2,] 113.69098+356.0602i -68.62964+200.0851i

# the following exmpales throw an error
if (FALSE) { # \dontrun{
rbind(a1, a2)   # the number of columns does not coincide 
cbind(pseries(a1), zvalues(a2)) # there is no automatic coercion
                                 # from pseries to zvalues
} # }
```
