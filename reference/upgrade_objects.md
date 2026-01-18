# Upgrade Objects to Common Class

Used in [`cbind`](https://rdrr.io/r/base/cbind.html),
[`rbind`](https://rdrr.io/r/base/cbind.html), `ratm_mult`, and
[`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
to make operations on different subclasses of the (super-) class `ratm`
possible.

## Usage

``` r
upgrade_objects(force = TRUE, ...)
```

## Arguments

- force:

  Default set to TRUE. If FALSE, then objects are also coerced when only
  one object is supplied. Option FALSE is used in bind methods

- ...:

  Arbitrary number of objects with superclass `ratm`

## Value

Upgraded inputs

## Details

Constant matrices are always upgraded to
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
and
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
and
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
are always upgraded to
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
if two `ratm` objects are involved. Therefore, most operations are
performed in the state space setting.

Depending on the highest "grade" of any object involved in the
operation, the lower ranked object is upgraded to the highest one. The
following ranks are assigned to the `ratm` objects:

- 1:
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)

- 2:
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)

- 3:
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)

- 4:
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)

- 5:
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)

- 6:
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)

- 7:
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)

Laurent polynomials
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
have a special place because they are detached from/unconnected to
elements
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
[`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md).

Operations with
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
are only possible with
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
objects (and those which can be coerced to
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)).

Operations which are not defined are

- [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  and
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md):
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  cannot be coerced (easily) to
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  object

- Anything involving
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  except for
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects

## Examples

``` r
p = test_polm(degree = 2)
lp = test_lpolm(degree_max = 1, degree_min = -1)
l = test_lmfd()
r = test_lmfd()
ss = test_stsp(s = 1)
ps = pseries(ss)
zv = zvalues(ss)

# Only one argument
rationalmatrices:::upgrade_objects(force = TRUE, l)
#> [[1]]
#> statespace realization [1,1] with s = 1 states
#>           s[1]       u[1]
#> s[1] 0.8016991 -0.4665939
#> x[1] 1.0000000  0.4609451
#> 
rationalmatrices:::upgrade_objects(force = FALSE, l)
#> [[1]]
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]   z^1 [,1]
#> [1,]        1 -0.8016991
#> right factor b(z):
#>       z^0 [,1]   z^1 [,1]
#> [1,] 0.4609451 -0.8361331
#> 

# Upgrade poly to rmfd
rationalmatrices:::upgrade_objects(force = TRUE, p, r)
#> [[1]]
#> statespace realization [1,1] with s = 2 states
#>      s[1] s[2] u[1]
#> s[1]    0    0    1
#> s[2]    1    0    0
#> x[1]  111  112  110
#> 
#> [[2]]
#> statespace realization [1,1] with s = 1 states
#>          s[1]      u[1]
#> s[1] 1.362329 -2.200286
#> x[1] 1.000000 -1.439654
#> 
rationalmatrices:::upgrade_objects(force = TRUE, l, ps)
#> [[1]]
#> ( 1 x 1 ) impulse response with maximum lag = 5 
#>      lag=0 [,1] lag=1 [,1] lag=2 [,1] lag=3 [,1] lag=4 [,1] lag=5 [,1]
#> [1,]  0.4609451 -0.4665939 -0.3740679 -0.2998899 -0.2404215 -0.1927457
#> 
#> [[2]]
#> ( 1 x 1 ) impulse response with maximum lag = 5 
#>      lag=0 [,1] lag=1 [,1] lag=2 [,1]  lag=3 [,1]  lag=4 [,1]  lag=5 [,1]
#> [1,]          1 -0.4762331  -0.210765 -0.09327765 -0.04128161 -0.01826988
#> 
```
