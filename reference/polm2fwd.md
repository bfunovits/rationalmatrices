# Transforms to Polynomial in Forward Shift

Transform
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
object to polynomial in forward shift (represented as
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object), i.e. transform \$\$a(z) = a_0 + a_1 z^1 + \cdots + a_p z^p\$\$
to \$\$a(z) = a_0 + a_1 z^{-1} + \cdots + a_p z^{-p}\$\$

## Usage

``` r
polm2fwd(polm_obj)
```

## Arguments

- polm_obj:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object

## Value

[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object

## Examples

``` r
(p = test_polm(degree = 3))
#> ( 1 x 1 ) matrix polynomial with degree <= 3 
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1]
#> [1,]      110      111      112      113
polm2fwd(p)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 0, and minimal degree >= -3
#>      z^-3 [,1] z^-2 [,1] z^-1 [,1] z^0 [,1]
#> [1,]       113       112       111      110
```
