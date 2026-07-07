# Coerce arrays to data frames

This helper function creates a `data.frame` from a given array.

## Usage

``` r
array2data.frame(x, rows = NULL, cols = NULL)
```

## Arguments

- x:

  Matrix or array.

- rows, cols:

  integer vectors. These two vectors define a partition of the
  "dimensions" (1,...,n), where n is the number of dimensions of x (i.e.
  length(dim(x))). If either of the two is missing, then the complement
  is used. At least one of the arguments "rows" and "cols" has to be
  given.

## Value

data.frame

## See also

The helper function `array2data.frame` is used internally for the
[`plot methods`](https://bfunovits.github.io/rationalmatrices/reference/plot.md).
The function
[`bmatrix`](https://bfunovits.github.io/rationalmatrices/reference/bmatrix.md)
(which coerces arrays to matrices) is a simplified version of
`array2data.frame`.

## Examples

``` r
# test array
x = test_array(dim = c(2,3,2), dimnames = TRUE)
array2data.frame(x, cols = c(1,3,2))
#>   A=1.C=1.B=1 A=2.C=1.B=1 A=1.C=2.B=1 A=2.C=2.B=1 A=1.C=1.B=2 A=2.C=1.B=2
#> 1         111         211         112         212         121         221
#>   A=1.C=2.B=2 A=2.C=2.B=2 A=1.C=1.B=3 A=2.C=1.B=3 A=1.C=2.B=3 A=2.C=2.B=3
#> 1         122         222         131         231         132         232
array2data.frame(x, rows = 1)
#>     A B=1.C=1 B=2.C=1 B=3.C=1 B=1.C=2 B=2.C=2 B=3.C=2
#> 1 A=1     111     121     131     112     122     132
#> 2 A=2     211     221     231     212     222     232
array2data.frame(x, rows = 2:1)
#>     B   A C=1 C=2
#> 1 B=1 A=1 111 112
#> 2 B=2 A=1 121 122
#> 3 B=3 A=1 131 132
#> 4 B=1 A=2 211 212
#> 5 B=2 A=2 221 222
#> 6 B=3 A=2 231 232
array2data.frame(x, rows = c(2,1,3))
#>      B   A   C value
#> 1  B=1 A=1 C=1   111
#> 2  B=2 A=1 C=1   121
#> 3  B=3 A=1 C=1   131
#> 4  B=1 A=2 C=1   211
#> 5  B=2 A=2 C=1   221
#> 6  B=3 A=2 C=1   231
#> 7  B=1 A=1 C=2   112
#> 8  B=2 A=1 C=2   122
#> 9  B=3 A=1 C=2   132
#> 10 B=1 A=2 C=2   212
#> 11 B=2 A=2 C=2   222
#> 12 B=3 A=2 C=2   232

# consider a pseudo socio economic data set
x = test_array(dim = c(2,4,5), random = TRUE)
dimnames(x) = list(sex=c('female','male'),
                   education = c('none','primary','high','university'),
                   age = c('<20','30-40','40-50','50-60','>60'))
array2data.frame(x, cols = 1)
#>     education   age      female        male
#> 1        none   <20 -0.85520250 -0.28689522
#> 2     primary   <20  0.89496163  0.06730444
#> 3        high   <20 -0.16267634 -0.82731017
#> 4  university   <20  1.87650562  0.76644020
#> 5        none 30-40  0.97995670  1.32178099
#> 6     primary 30-40 -1.11971083  0.51459982
#> 7        high 30-40 -1.50909984  1.53274148
#> 8  university 30-40  0.42914737  0.12210341
#> 9        none 40-50 -1.13801240 -0.55801513
#> 10    primary 40-50  1.05253854  0.67768364
#> 11       high 40-50  0.03849955 -0.35638119
#> 12 university 40-50  0.78284410  0.80441162
#> 13       none 50-60 -1.90006082  0.93578429
#> 14    primary 50-60 -0.30905150  0.26306668
#> 15       high 50-60 -1.79059186 -0.78825884
#> 16 university 50-60 -1.13302167  0.36365257
#> 17       none   >60 -0.28588791  0.51766913
#> 18    primary   >60 -0.10290867 -0.97406959
#> 19       high   >60  1.27067230  0.96086479
#> 20 university   >60  0.76872137  1.03593077
array2data.frame(x, cols = c(1,2))
#>     age female.none  male.none female.primary male.primary female.high
#> 1   <20  -0.8552025 -0.2868952      0.8949616   0.06730444 -0.16267634
#> 2 30-40   0.9799567  1.3217810     -1.1197108   0.51459982 -1.50909984
#> 3 40-50  -1.1380124 -0.5580151      1.0525385   0.67768364  0.03849955
#> 4 50-60  -1.9000608  0.9357843     -0.3090515   0.26306668 -1.79059186
#> 5   >60  -0.2858879  0.5176691     -0.1029087  -0.97406959  1.27067230
#>    male.high female.university male.university
#> 1 -0.8273102         1.8765056       0.7664402
#> 2  1.5327415         0.4291474       0.1221034
#> 3 -0.3563812         0.7828441       0.8044116
#> 4 -0.7882588        -1.1330217       0.3636526
#> 5  0.9608648         0.7687214       1.0359308

# convert an impulse response function into a data.frame
# generate a random statespace model with m=3 outputs, n=2 inputs and s=4 states.
model = stsp(A = matrix(stats::rnorm(16), nrow = 4, ncol = 4), 
             B = matrix(stats::rnorm(8),  nrow = 4, ncol = 2),
             C = matrix(stats::rnorm(12), nrow = 3, ncol = 4), 
             D = matrix(stats::rnorm(6), nrow = 3, ncol = 2))
k = unclass(pseries(model, lag.max = 25))
dimnames(k) = list(y = paste('y[',1:3,']',sep=''), 
                   x = paste('x[',1:2,']',sep=''),
                   lag = 0:25)
head(array2data.frame(k, rows = c(1,3), cols=2))
#>      y lag       x[1]       x[2]
#> 1 y[1]   0  0.7863626  0.0751059
#> 2 y[2]   0 -1.2705131  0.5585144
#> 3 y[3]   0  0.5421415  0.4154064
#> 4 y[1]   1 -1.0742136 -0.2137160
#> 5 y[2]   1 -0.9738554 -3.7657659
#> 6 y[3]   1 -3.6647290 -2.7221471
```
