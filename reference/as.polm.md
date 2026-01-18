# Coerece to polynomial object

Coerece to polynomial object

## Usage

``` r
as.polm(obj, ...)

# S3 method for class 'lpolm'
as.polm(obj, ...)
```

## Arguments

- obj:

  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  object

- ...:

  other arguments

## Value

[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
object

## Examples

``` r
lp = lpolm(1:3, min_deg = 1) 
as.polm(lp)
#> ( 1 x 1 ) matrix polynomial with degree <= 3 
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1]
#> [1,]        0        1        2        3
```
