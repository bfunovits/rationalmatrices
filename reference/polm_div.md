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
#>       z^0 [,1]      [,2]  z^1 [,1]      [,2]  z^2 [,1]      [,2]
#> [1,]  36.12466  44.29911 13.449209 14.544444 -5.466176 -6.566854
#> [2,] -65.99837 -72.09624  5.327070  7.312464  6.349676  6.154500
#> [3,]  71.40729  81.53618  7.264565  6.581119 -8.078047 -8.931964
#> 
#> $rem
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>        z^0 [,1]      [,2]  z^1 [,1]        [,2]
#> [1,]   60.46469 -16.82513  46.56593  -0.9500461
#> [2,] -100.45128  18.90242 -19.49479 -24.3284306
#> [3,]  111.57581 -21.92424  46.10258  14.3132742
#> 

all.equal(a, out$qu %r% b + out$rem)
#> [1] TRUE
```
