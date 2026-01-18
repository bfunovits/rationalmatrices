# Column End Matrix of a Polynomial Matrix

The *column end matrix* of an \\(m,n)\\-dimensional polynomial matrix
\\a(z)=a_0 + a_1 z + \cdots + a_p z^p\\ is defined as follows. Suppose
that the maximum degree of the elements in the \\i\\-th column is
\\p_i\\. Then the column end matrix is the \\(m,n)\\ matrix with
\\i\\-th column equal to the \\i\\-th column of the coefficient matrix
\\a\_{p_i}\\. If a column of \\a(z)\\ is zero, then the elements of the
corresponding column of the column end matrix are set to `NA`'s.

## Usage

``` r
col_end_matrix(x)
```

## Arguments

- x:

  A polynomial matrix, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

## Value

The column end matrix.

## Examples

``` r
x = polm(array(c(0,1,1,0,
                 0,0,1,0,
                 0,0,0,1,
                 0,0,0,0), dim = c(2,2,4)))
x
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]        0     1        0     1        0     0        0     0
#> [2,]        1     0        0     0        0     1        0     0
degree(x)
#>      [,1] [,2]
#> [1,]   -1    1
#> [2,]    0    2
degree(x, 'columns')
#> [1] 0 2
col_end_matrix(x)
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    1    1
```
