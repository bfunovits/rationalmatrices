# Coerce to Laurent polynom object

The attribute `min_deg` is set to zero for the given function input.

## Usage

``` r
as.lpolm(obj, ...)

# S3 method for class 'polm'
as.lpolm(obj, ...)
```

## Arguments

- obj:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object

- ...:

  other arguments

## Value

[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object

## Examples

``` r
p = test_polm(degree = 2)
as.lpolm(p)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 2, and minimal degree >= 0
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]      110      111      112
```
