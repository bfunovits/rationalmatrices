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
#>      z^0 [,1]  [,2]   z^1 [,1]        [,2]
#> [1,]        1     0 -0.2648497 -0.04616607
#> [2,]        0     1  0.3708954 -0.23294641
#> right factor b(z):
#>         z^0 [,1]       [,2]   z^1 [,1]      [,2]
#> [1,] -1.05570358  0.3521215 -0.2102821 0.7562880
#> [2,] -0.05883459 -0.4889416 -0.1687958 0.5326999
zvalues(ll)
#> ( 2 x 2 ) frequency response
#>          z[1] [,1]          [,2]             z[2] [,1]                  [,2]
#> [1,] -1.6894143+0i  1.4667756+0i -1.0921742+0.5252345i  0.3805187-0.8860062i
#> [2,]  0.5201274+0i -0.6521864+0i -0.2437643-0.2488275i -0.1224071-0.2626423i
#>                   z[3] [,1]                  [,2]              z[4] [,1]
#> [1,] -0.6995270+0.19406056i -0.2350832-0.3088891i -0.6995270-0.19406056i
#> [2,] -0.1444185+0.02079844i -0.8102761-0.2912206i -0.1444185-0.02079844i
#>                       [,2]             z[5] [,1]                  [,2]
#> [1,] -0.2350832+0.3088891i -1.0921742-0.5252345i  0.3805187+0.8860062i
#> [2,] -0.8102761+0.2912206i -0.2437643+0.2488275i -0.1224071+0.2626423i
```
