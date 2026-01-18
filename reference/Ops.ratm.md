# Arithmetic Ops Group Methods for Rational Matrices

Implements the following basic arithmetic operations on rational
matrices

- unary operators `'+a'` and `'-a'`.

- the power (`'a^k'`) operator for *square, non empty* rational matrices
  `a` and integer powers `k`.

- elementwise multiplication (`'a * b'`)

- addition and substraction (`'a + b'` and `'a - b'`)

- elementwise polynomal division (`'a %/% b'`),

- elementwise polynomial remainder (`'a %% b'`),

## Usage

``` r
# S3 method for class 'ratm'
Ops(e1, e2)
```

## Arguments

- e1, e2:

  At least one of `e1, e2` must be an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md).

## Value

Rational matrix object.

## Details

The unitary operators `'+a'` and `'-a'` are implemented for all classes.

The power (`'a^k'`) operator is only implemented *square*, rational
matrices and integer powers `k`.

- `a^0` returns the identity matrix (represented by an object of the
  same class as the input argument `a`).

- `a^1` simply returns the input arguments `a`.

- The case \\k\>1\\ is implemented for all classes. However, `lmfd` and
  `rmfd` objects are first coerced to `stsp` objects.

- The case \\k\<0\\ is implemented for all classes, except for Laurent
  polynomials (`lpolm` objects). However, the matrix must be non empty
  and `polm`, `lmfd` and `rmfd` objects are first coerced to `stsp`
  objects.

For the binary operators \`a + b\`, \`a - b\` and \`a \* b\` the two
arguments are first coerced to a common class according to the following
scheme

|                |     |         |                                                  |
|----------------|-----|---------|--------------------------------------------------|
| matrix         | -\> | polm    | (coerce non rational matrices to `polm` objects) |
| lmfd           | -\> | stsp    | (coerce `lmfd` objects to `stsp` objects)        |
| rmfd           | -\> | stsp    | (coerce `rmfd` objects to `stsp` objects)        |
| polm o stsp    | -\> | stsp    |                                                  |
| polm o pseries | -\> | pseries |                                                  |
| polm o zvalues | -\> | zvalues |                                                  |
| stsp o pseries | -\> | pseries |                                                  |
| stsp o zvalues | -\> | zvalues |                                                  |

Note that `pseries` objects cannot be (easily) coerced to `zvalues`
objects, so these binary operations throw an error, if one tries to
combine a `pseries` with a `zvalues` object.

If two `pseries` objects are combined then they are truncated to the
minimum of the respective number of “lags”. Two `zvalues` objects are
only combined if the “z” values are identical. Otherwise an error is
thrown.

Of course the two arguments have to be of compatible dimension. If one
of the arguments is a scalar (\\(1,1)\\ matrix) then this argument is
"expanded" to a compatible matrix with identical entries.

Note that the computed statespace realizations are often non minimal!
(This remark also applies for other operations on statespace
realizations.)

For the matrix multiplication, see
[`%r%`](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md).

The elementwise polynomal division (`'a %/% b'`) and the elementwise
polynomial remainder (`'a %% b'`) are of course only implemented for
polynomial matrices (`polm` objects) or (objects which may be coerced to
`polm` objects).

The above remark on scalar arguments also applies for these operations.

## Examples

``` r
# Multiplication (and division) of a scalar (from left and right)
a = test_polm(dim = c(2,2), degree = 3)
a
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]      110   120      111   121      112   122      113   123
#> [2,]      210   220      211   221      212   222      213   223
-a
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]     -110  -120     -111  -121     -112  -122     -113  -123
#> [2,]     -210  -220     -211  -221     -212  -222     -213  -223
a * (-1)
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]     -110  -120     -111  -121     -112  -122     -113  -123
#> [2,]     -210  -220     -211  -221     -212  -222     -213  -223
a %/% 0.5
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]      220   240      222   242      224   244      226   246
#> [2,]      420   440      422   442      424   444      426   446

# Addition
2 + a        # 2 is coerced to a constant matrix polynomial
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]      112   122      111   121      112   122      113   123
#> [2,]      212   222      211   221      212   222      213   223

# Elementwise remainder
a %% 0.5
#> ( 2 x 2 ) matrix polynomial with degree <= -1 
0.5 %% a
#> ( 2 x 2 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]  [,2]
#> [1,]      0.5   0.5
#> [2,]      0.5   0.5

# Elementwise division and multiplication with scalar polm
z = polm(c(0,1))
a %/% z
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2]
#> [1,]      111   121      112   122      113   123
#> [2,]      211   221      212   222      213   223
z * a
#> ( 2 x 2 ) matrix polynomial with degree <= 4 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2] z^4 [,1]  [,2]
#> [1,]        0     0      110   120      111   121      112   122      113   123
#> [2,]        0     0      210   220      211   221      212   222      213   223
a * z^2
#> ( 2 x 2 ) matrix polynomial with degree <= 5 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2] z^4 [,1]  [,2]
#> [1,]        0     0        0     0      110   120      111   121      112   122
#> [2,]        0     0        0     0      210   220      211   221      212   222
#>      z^5 [,1]  [,2]
#> [1,]      113   123
#> [2,]      213   223

# (Non-negative integer) power of univariate polynomial 
# (useful for generating a polynomial matrix)
z^3
#> ( 1 x 1 ) matrix polynomial with degree <= 3 
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1]
#> [1,]        0        0        0        1
matrix(0, nrow = 2) + z * matrix(1, nrow = 2) + z^2 * matrix(2, nrow = 2)
#> ( 2 x 1 ) matrix polynomial with degree <= 2 
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]        0        1        2
#> [2,]        0        1        2

# (Non-negative integer) power of quadratic polynomial matrices
a^0
#> ( 2 x 2 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
a^1
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]      110   120      111   121      112   122      113   123
#> [2,]      210   220      211   221      212   222      213   223
a^2
#> ( 2 x 2 ) matrix polynomial with degree <= 6 
#>      z^0 [,1]  [,2] z^1 [,1]   [,2] z^2 [,1]   [,2] z^3 [,1]   [,2] z^4 [,1]
#> [1,]    37300 39600    75150  79770   113552 120512   152508 161828   115220
#> [2,]    69300 73600   139350 147970   210152 223112   281708 299028   212420
#>        [,2] z^5 [,1]   [,2] z^6 [,1]  [,2]
#> [1,] 122240    77374  82074    38968 41328
#> [2,] 225440   142374 151074    71568 75928

# inverse of square polynomial matrix
a = test_polm(dim = c(2,2), degree = 1, random = TRUE, digits = 2)
a^(-1)
#> statespace realization [2,2] with s = 2 states
#>          s[1]      s[2]     u[1]     u[2]
#> s[1] 0.350000 -4.162500 1.875000 8.750000
#> s[2] 1.621429 -1.844643 2.053571 4.821429
#> x[1] 0.350000 -4.162500 1.875000 8.750000
#> x[2] 1.621429 -1.844643 2.053571 4.821429
print(pseries(a %r% a^(-1)), digits = 3)
#> ( 2 x 2 ) impulse response with maximum lag = 5 
#>      lag=0 [,1]  [,2] lag=1 [,1]  [,2] lag=2 [,1]  [,2] lag=3 [,1]  [,2]
#> [1,]          1     0          0     0          0     0          0     0
#> [2,]          0     1          0     0          0     0          0     0
#>      lag=4 [,1]  [,2] lag=5 [,1]  [,2]
#> [1,]          0     0          0     0
#> [2,]          0     0          0     0

# elementwise multiplication 
a = test_polm(dim = c(2,3), degree = 1, random = TRUE, digits = 2)
b = test_polm(dim = c(2,3), degree = 1, random = TRUE, digits = 2)
b * a
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]    [,2]   [,3] z^1 [,1]    [,2]    [,3] z^2 [,1]    [,2]    [,3]
#> [1,]  -0.5538 -0.3192 0.9095  -2.9073 -0.6743 -0.2933  -0.5895  0.1935 -0.0058
#> [2,]  -0.2226 -0.2875 0.1632  -1.1315 -0.3415  1.8636  -1.3066 -0.0946  0.5340
all.equal(zvalues(a * b), zvalues(a) * zvalues(b))
#> [1] TRUE
```
