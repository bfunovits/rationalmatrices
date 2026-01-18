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
#>       z^0 [,1]       [,2]    z^1 [,1]      [,2]   z^2 [,1]       [,2]
#> [1,] -1.947685  -7.115806  0.03724353  1.660189 -0.2301326 -0.3610648
#> [2,] -4.915987 -14.650755  1.22579629  3.471123 -0.5963417 -1.3758817
#> [3,]  6.534776  21.122758 -2.42079170 -5.812314  0.8515300  2.1146335
#> 
#> $rem
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>       z^0 [,1]      [,2]   z^1 [,1]      [,2]
#> [1,]  2.277420 -11.61503   5.849339  11.02783
#> [2,]  1.935910 -27.33473  15.568595  22.96784
#> [3,] -2.932782  38.89806 -20.845248 -35.56956
#> 

all.equal(a, out$qu %r% b + out$rem)
#> [1] TRUE
```
