# QL and LQ Decomposition

Returns the QL and LQ decomposition of a matrix with non-negative
"diagonal" elements in the QL and LQ decompositions. Only works if no
column pivoting occurs in the QR decomposition
[`qr`](https://rdrr.io/r/base/qr.html).

## Usage

``` r
ql_decomposition(x, ...)

lq_decomposition(x, ...)
```

## Arguments

- x:

  Matrix.

- ...:

  Other arguments for [`qr`](https://rdrr.io/r/base/qr.html)

## Value

List with two elements:

- q:

  Semi-orthogonal matrix

- l:

  Lower triangular matrix with non-negative diagonal elements.

## Note

Shouldn't export this function when publishing package!

## Implementation of the QL decomposition using the QR decomposition

The base function [qr](https://rdrr.io/r/base/qr.html) is used in the
following way. (We need to assume that there is no pivoting since
otherwise the function throws an error). First, we reorder the columns
of the input matrix `x` from last to first QR-decompose this matrix.
Next, we reorder the columns of `Q` and the rows of `R` in the same way.
The original matrix `x` is now in `QL` decomposition.

## Examples

``` r
set.seed(1803)

# Tall matrix
x = matrix(stats::rnorm(5*3), 5, 3)
out = ql_decomposition(x)
all.equal(x, out$q %*% out$l)
#> [1] TRUE
all.equal(diag(ncol(out$q)), t(out$q) %*% out$q)
#> [1] TRUE
out$l
#>            [,1]     [,2]     [,3]
#> [1,]  2.3430442 0.000000 0.000000
#> [2,] -0.4813528 2.504091 0.000000
#> [3,]  0.2570508 1.183947 2.350825
out = lq_decomposition(x)
all.equal(x, out$l %*% out$q)
#> [1] TRUE
all.equal(diag(nrow(out$q)), out$q %*% t(out$q))
#> [1] TRUE
out$l
#>            [,1]       [,2]       [,3]
#> [1,]  2.4318488  0.0000000  0.0000000
#> [2,]  0.1373582  1.0693912  0.0000000
#> [3,] -0.9217318 -1.7489450  0.7996942
#> [4,]  0.6810998 -0.7241543 -2.0071689
#> [5,] -0.8357022 -1.1436741 -0.5818361

# Wide matrix
x = matrix(stats::rnorm(5*3), 3, 5)
out = ql_decomposition(x)
all.equal(x, out$q %*% out$l)
#> [1] TRUE
all.equal(diag(ncol(out$q)), t(out$q) %*% out$q)
#> [1] TRUE
out$l
#>            [,1]       [,2]       [,3]     [,4]     [,5]
#> [1,]  1.0854922 1.16839176  0.6295948 0.000000 0.000000
#> [2,]  1.0702433 1.45849417  0.1366193 1.302877 0.000000
#> [3,] -0.8729544 0.08863472 -0.1230865 1.160753 2.179503
out = lq_decomposition(x)
all.equal(x, out$l %*% out$q)
#> [1] TRUE
all.equal(diag(nrow(out$q)), out$q %*% t(out$q))
#> [1] TRUE
out$l
#>            [,1]      [,2]     [,3]
#> [1,]  1.4032852 0.0000000 0.000000
#> [2,] -0.2164448 2.6480266 0.000000
#> [3,]  1.8590980 0.2003933 1.512196

# Square matrix
x = matrix(stats::rnorm(4*4), 4, 4)
out = ql_decomposition(x)
all.equal(x, out$q %*% out$l)
#> [1] TRUE
all.equal(diag(ncol(out$q)), t(out$q) %*% out$q)
#> [1] TRUE
out$l
#>            [,1]       [,2]       [,3]     [,4]
#> [1,]  0.7835124  0.0000000  0.0000000 0.000000
#> [2,]  1.6198819  0.4521451  0.0000000 0.000000
#> [3,] -0.5476521 -0.1655997  2.3140488 0.000000
#> [4,]  0.4083269 -1.1031379 -0.1874477 2.482612
out = lq_decomposition(x)
all.equal(x, out$l %*% out$q)
#> [1] TRUE
all.equal(diag(nrow(out$q)), out$q %*% t(out$q))
#> [1] TRUE
out$l
#>           [,1]       [,2]      [,3]      [,4]
#> [1,] 2.1763716  0.0000000  0.000000 0.0000000
#> [2,] 0.2516972  2.0028011  0.000000 0.0000000
#> [3,] 0.7376794 -1.5175622  1.353042 0.0000000
#> [4,] 1.0090988 -0.5026001 -1.351903 0.3450825

# reset seed
set.seed(NULL)
```
