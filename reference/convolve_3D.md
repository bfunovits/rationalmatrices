# Convolution of Matrix-valued Sequences

![Internal function](figures/internal.svg)

Compute the convolution of two (matrix valued) sequences \\a_0, a_1,
\ldots\\ and \\b_0, b_1, \ldots\\, i.e. \$\$c_k = \sum\_{j=0}^{k} a_j
b\_{k-j}\$\$

## Usage

``` r
convolve_3D(a, b, truncate = FALSE)
```

## Arguments

- a:

  3D array with dimension `(m,n,p+1)`

- b:

  3D array with dimension `(n,o,q+1)`

- truncate:

  (boolean) if `TRUE` the output sequence has length \\\min(p,q)+1\\,
  otherwise the sequence has length \\p+q+1\\.

## Value

3D array of dimension `(m,o,r+1)`, where `r=p+q` if `truncate==FALSE`
and `r = min(p,q)` otherwise.

## See also

This helper function is used for the multiplication of `polm`, `lpolm`
and `pseries` objects, see
[`%r%`](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md).

## Examples

``` r
a = test_array(dim = c(3,2,2), random = TRUE)
b = test_array(dim = c(2,1,3), random = TRUE)
convolve_3D(a,b)
#> , , 1
#> 
#>            [,1]
#> [1,] -0.1847338
#> [2,] -0.5227425
#> [3,] -0.2158141
#> 
#> , , 2
#> 
#>            [,1]
#> [1,] -0.6699038
#> [2,] -3.9129920
#> [3,] -0.2122951
#> 
#> , , 3
#> 
#>           [,1]
#> [1,] -3.372669
#> [2,]  2.183086
#> [3,] -1.354943
#> 
#> , , 4
#> 
#>            [,1]
#> [1,] -0.5285895
#> [2,]  0.3638544
#> [3,] -0.4597577
#> 
convolve_3D(a,b,TRUE)
#> , , 1
#> 
#>            [,1]
#> [1,] -0.1847338
#> [2,] -0.5227425
#> [3,] -0.2158141
#> 
#> , , 2
#> 
#>            [,1]
#> [1,] -0.6699038
#> [2,] -3.9129920
#> [3,] -0.2122951
#> 
```
