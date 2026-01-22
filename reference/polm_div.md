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
#>        z^0 [,1]      [,2]  z^1 [,1]      [,2]  z^2 [,1]      [,2]
#> [1,]   68.63017  133.9147  10.11417  19.19290  2.182392  3.369295
#> [2,]  302.39270  590.3264  43.09670  76.68699  4.310962 13.347594
#> [3,] -214.09961 -410.3890 -27.34128 -56.95725 -3.011508 -6.213239
#> 
#> $rem
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>        z^0 [,1]       [,2]   z^1 [,1]      [,2]
#> [1,]  -40.74383  -54.18968  -78.15372  240.8611
#> [2,] -182.98575 -241.99021 -348.40940 1067.6433
#> [3,]  127.80225  169.94282  245.16864 -743.3890
#> 

all.equal(a, out$qu %r% b + out$rem)
#> [1] TRUE
```
