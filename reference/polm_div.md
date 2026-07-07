# Division Algorithm for Polynomial Matrices

For given polynomial matrices \\a(z), b(z)\\ compute two matrices
\\c(z), d(z)\\ such that \$\$a(z) = c(z) b(z) + d(z)\$\$ where the
degree of \\d(z)\\ is smaller than the degree of \\b(z)\\. The matrix
\\b(z)\\ must be square with a non singular leading coefficient matrix!
The matrices must be compatible, i.e. the number of columns of \\a(z)\\
must equal the number of rows (and columns) of \\b(z)\\.

## Usage

``` r
polm_div(a, b)
```

## Arguments

- a, b:

  Two compatible polynomial matrices.

## Value

List with two slots

- qucontains the polynomial \\c(z)\\

- remcontains the polynomial \\d(z)\\.

## Examples

``` r
a = test_polm(dim = c(3,2), degree = 4, random = TRUE)
b = test_polm(dim = c(2,2), degree = 2, random = TRUE)

(out = polm_div(a, b))
#> $qu
#> ( 3 x 2 ) matrix polynomial with degree <= 2 
#>        z^0 [,1]       [,2]    z^1 [,1]        [,2]  z^2 [,1]       [,2]
#> [1,] -1.1407040  1.3689253 -0.05007305 -0.98069292 0.1292603  0.9231022
#> [2,] -0.2694931 -0.2119719  0.59082032 -0.04523667 0.2019229  0.7190641
#> [3,]  4.2526376  1.2255002 -0.72896475  2.62352494 0.2758249 -1.1449855
#> 
#> $rem
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>        z^0 [,1]       [,2]   z^1 [,1]      [,2]
#> [1,] -0.5390637  2.6952921  4.4251373  1.432490
#> [2,]  1.3147424 -0.1589452  0.7688883 -2.105304
#> [3,]  3.0560400 -3.3818161 -6.3485225  2.066395
#> 

all.equal(a, out$qu %r% b + out$rem)
#> [1] TRUE
```
