# Frequency Response Function

This function evaluates a rational matrix at given (complex) arguments.
If the rational matrix corresponds to a *linear, dynamic filter*, then
this rational function is also called *transfer function* of the filter
and the values of the function on the complex unit circle are called
*frequency response* of the filter.

## Usage

``` r
zvalues(obj, z, f, n.f, sort.frequencies, ...)

# Default S3 method
zvalues(obj, z = NULL, f = NULL, n.f = NULL, sort.frequencies = FALSE, ...)

# S3 method for class 'polm'
zvalues(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE, ...)

# S3 method for class 'lpolm'
zvalues(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE, ...)

# S3 method for class 'lmfd'
zvalues(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE, ...)

# S3 method for class 'rmfd'
zvalues(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE, ...)

# S3 method for class 'stsp'
zvalues(obj, z = NULL, f = NULL, n.f = 5, sort.frequencies = FALSE, ...)
```

## Arguments

- obj:

  (rational) matrix object, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object. or an object which may be coerced to a polynomial matrix with
  `polm(obj)`. The default `S3` method acts as a constructor for zvalues
  objects: and array object and `z` are supplied, corresponding to the
  evaluated object at the values of `z`. (Numeric or complex) vectors or
  matrices are coerced arrays. If `f` is supplied instead of `z`, then
  `z` is considered to be `exp(complex(imaginary = (-2*pi)*f))`.

- z:

  (numeric or complex) vector, points at which to evaluate the rational
  matrix. Theses values are stored as attribute and can be accessed
  through `zvalues_obj$z`.

- f:

  (numeric) vector of frequencies. If `z = NULL` then `z` is set to
  `z = exp(complex(imaginary = (-2*pi)*f))`. If `z` is non `NULL`, then
  `f` is set to `f = -Arg(z)/(2*pi)`. Theses values can be accessed
  through `zvalues_obj$f`.

- n.f:

  (integer) number of frequencies. If `z = f = NULL` then a grid of
  frequencies `f = (0:(n.f-1))/n.f` is used (and `z` is generated as
  explained above).

- sort.frequencies:

  boolean, sort frequencies

- ...:

  optional additional parameters

## Value

Object of type `zvalues`

## Details

A `zvalues` object is simply a (complex valued) 3-dimensional array of
dimension \\(m,n,n.f)\\ with attribute `z` and class attribute
`c('zvalues','ratm')`. The dimensions \\(m,n,n.f)\\ may also be zero,
where `n.f` is the length of `z`.

`z` and `f` can be accessed through `zvalues_obj$z` and `zvalues_obj$z`,
respectively.

If you only need the value of the matrix evaluated at one point, then
use
[`zvalue`](https://bfunovits.github.io/rationalmatrices/reference/zvalue.md)`(obj, z)`.

## Examples

``` r
(vv = zvalues(1:3, z = 4:6))
#> ( 1 x 1 ) frequency response
#>      z[1] [,1] z[2] [,1] z[3] [,1]
#> [1,]         1         2         3
(mm = zvalues(matrix(1:4, 2, 2), z = 1))
#> ( 2 x 2 ) frequency response
#>      z[1] [,1]  [,2]
#> [1,]         1     3
#> [2,]         2     4
(aa = zvalues(array(1:8, dim = c(2,2,2)), f = c(0.2, 0.8)))
#> ( 2 x 2 ) frequency response
#>      z[1] [,1]  [,2] z[2] [,1]  [,2]
#> [1,]         1     3         5     7
#> [2,]         2     4         6     8
(ll = test_lmfd(dim = c(2,2), bpoles = 1, bzeroes = 1))
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2]    z^1 [,1]        [,2]
#> [1,]        1     0 -0.01873117 -0.08245998
#> [2,]        0     1 -0.07300167 -0.46444409
#> right factor b(z):
#>        z^0 [,1]       [,2]  z^1 [,1]       [,2]
#> [1,]  0.3345009 -0.4157322 0.2438196  0.8044345
#> [2,] -0.2158452 -0.6666081 0.2679840 -0.2180263
zvalues(ll)
#> ( 2 x 2 ) frequency response
#>         z[1] [,1]          [,2]             z[2] [,1]                  [,2]
#> [1,] 0.6044650+0i  0.2602959+0i  0.3858013-0.2242859i -0.1526344-0.7028060i
#> [2,] 0.1797492+0i -1.6163249+0i -0.2654003-0.1978830i -0.6289792+0.5603479i
#>                 z[3] [,1]                  [,2]            z[4] [,1]
#> [1,]  0.1535565-0.123745i -1.0270949-0.4542050i  0.1535565+0.123745i
#> [2,] -0.3343848-0.047620i -0.2864177+0.2015201i -0.3343848+0.047620i
#>                       [,2]             z[5] [,1]                  [,2]
#> [1,] -1.0270949+0.4542050i  0.3858013+0.2242859i -0.1526344+0.7028060i
#> [2,] -0.2864177-0.2015201i -0.2654003+0.1978830i -0.6289792-0.5603479i
```
