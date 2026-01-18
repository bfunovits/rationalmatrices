# Ho-Kalman Realization Algorithm

This helper function implements the Ho-Kalman algorithm.

## Usage

``` r
pseries2stsp(
  obj,
  method = c("balanced", "echelon"),
  Hsize = NULL,
  s = NULL,
  nu = NULL,
  tol = sqrt(.Machine$double.eps),
  Wrow = NULL,
  Wcol = NULL
)
```

## Arguments

- obj:

  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  object or 3-D array with dimension \\(m,n,l+1)\\.

- method:

  Character string, which determines the method and the
  "parametrization" type of the state space model. See below for more
  details.

- Hsize:

  integer vector `c(f,p)`, number of block rows and block columns of the
  Hankel matrix which is used to construct the statespace realization.
  If NULL a default choice is made.

- s:

  desired state dimension. Only used for `method = "balanced"`. Note
  however, if \\s\\ is larger than the rank of the Hankel matrix, then
  the procedure will break down. If `s` is missing, then the state
  dimension is determined from the singular values of the Hankel matrix.
  To be precise the state dimension is chosen as the number of singular
  values which are greater than or equal to `tol` times the maximum
  singular value.

- nu:

  Kronecker indices. Only used for `method = "echelon"`. If missing,
  then `nu` is computed with a QR decomposition of the transpose of the
  Hankel matrix of the impulse response coefficients.

- tol:

  tolerance parameter used for the QR decomposition or the SVD
  decomposition of the Hankel matrix \\H\\ of the impulse response
  coefficients.

- Wrow, Wcol:

  weighting matrices (default is no weighting, i.e. identity matrices).
  These weighting matrices are only used for `method="balanced"`, where
  the SVD of the weighted Hankel matrix `Wrow %*% H %*% t(Wcol)` is
  computed.

## Value

List with slots

- Xs:

  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object, the rational matrix in statespace form

- Hsv:

  Singular values of the Hankel matrix for `method='balanced'` and
  `NULL` else.

- nu:

  Kronecker indices for `method='echelon'` and `NULL` else.

## Details

The procedure(s) may be used for model reduction (with some care).

There are a number of restrictions on the number of lags \\l\\ of the
impulse response, the number of block rows (\\f\\), block columns
(\\p\\) of the Hankel matrix and the Kronecker indices \\\nu_i\\. We
require that: \\p\>0\\, \\f\>1\\, \\l \geq f+p-1\\ and \\\nu_i \<f\\. If
these restrictions are not satisfied an error is thrown.

## Examples

``` r
# generate random rational matrix X(z) in statespace form
# make sure that the A matrix is stable
m = 3
n = 2
s = 7
A = matrix(rnorm(s*s), nrow = s, ncol = s)
A = A / (1.1 * max(abs(eigen(A, only.values = TRUE)$values)))
Xs = stsp(A, B = matrix(rnorm(s*n), nrow = s, ncol = n),
          C = matrix(rnorm(s*m), nrow = m, ncol = s),
          D = diag(1, nrow = m, ncol = n))
Xi = pseries(Xs, lag.max = 20)

out = pseries2stsp(Xi, method = 'balanced')
print(out)
#> $Xs
#> statespace realization [3,2] with s = 7 states
#>              s[1]         s[2]        s[3]         s[4]         s[5]
#> s[1] -0.923236142 -0.004059048 -0.18926223 -0.009955537  0.007992984
#> s[2] -0.167377808  0.165894886  0.12815816 -0.073323053 -0.040340847
#> s[3] -0.052770213  0.674825358 -0.06148503 -0.221446913 -0.139445970
#> s[4]  0.042550415  0.317309275 -0.18229272  0.695158655  0.376458820
#> s[5]  0.034320234  0.056547553 -0.37824989 -0.412841898  0.141875367
#> s[6]  0.022784855  0.116926181  0.18998354 -0.163395859 -0.003732914
#> s[7] -0.005555881 -0.012164435  0.04289856 -0.046149932  0.403363218
#> x[1] -1.455906816  0.603953058  1.49738378  0.114720965  0.099254675
#> x[2] -0.567885821 -0.800199631  0.25177463  0.556776459 -0.621187470
#> x[3]  0.444535363  1.121674640 -0.50482169  0.499017699 -0.385939590
#>              s[6]          s[7]        u[1]        u[2]
#> s[1] -0.001032311 -0.0040303712  0.65140060 -1.52907755
#> s[2] -0.190197502  0.1337757833 -1.69262697  0.77339130
#> s[3]  0.170254112  0.1050559752  0.75807008  0.87744489
#> s[4]  0.090875757  0.0986579737 -0.24278098 -0.42646482
#> s[5] -0.021500536 -0.1257277854 -0.19775104 -0.22137960
#> s[6]  0.102980559 -0.0007722748 -0.07446637 -0.31664263
#> s[7] -0.345775863  0.3068792650  0.05074257  0.05131137
#> x[1]  0.006923784 -0.1325445482  1.00000000  0.00000000
#> x[2]  0.096303165  0.1074984341  0.00000000  1.00000000
#> x[3] -0.210420361 -0.1900921386  0.00000000  0.00000000
#> 
#> $Hsv
#>  [1] 1.584987e+01 4.091085e+00 3.379543e+00 1.797861e+00 9.288931e-01
#>  [6] 3.429737e-01 2.298573e-01 1.142788e-15 8.663900e-16 7.177996e-16
#> [11] 6.067676e-16 4.436562e-16 3.768441e-16 3.215026e-16 2.851432e-16
#> [16] 2.558311e-16 2.030134e-16 1.727346e-16 1.447868e-16 1.237507e-16
#> 
#> $nu
#> NULL
#> 
# check impulse response
all.equal(pseries(out$Xs, lag.max = 20), Xi)
#> [1] TRUE

Xs1 = as.stsp(Xi)
all.equal(Xs1, out$Xs)
#> [1] TRUE

out = pseries2stsp(Xi, method = 'echelon')
print(out)
#> $Xs
#> statespace realization [3,2] with s = 7 states
#>              s[1]        s[2]        s[3]       s[4]       s[5]       s[6]
#> s[1]  0.000000000  0.00000000  0.00000000  1.0000000  0.0000000  0.0000000
#> s[2]  0.000000000  0.00000000  0.00000000  0.0000000  1.0000000  0.0000000
#> s[3]  0.000000000  0.00000000  0.00000000  0.0000000  0.0000000  1.0000000
#> s[4]  0.000000000  0.00000000  0.00000000  0.0000000  0.0000000  0.0000000
#> s[5] -0.098535478 -0.34250460  0.33115814 -0.8865510  2.0175720 -1.3090509
#> s[6] -0.105138543 -0.25605263 -0.07161845  0.1636435  0.2679175  0.6665081
#> s[7]  0.004956274  0.09768944 -0.65541641  1.6749946 -2.7121341  2.9281574
#> x[1]  1.000000000  0.00000000  0.00000000  0.0000000  0.0000000  0.0000000
#> x[2]  0.000000000  1.00000000  0.00000000  0.0000000  0.0000000  0.0000000
#> x[3]  0.000000000  0.00000000  1.00000000  0.0000000  0.0000000  0.0000000
#>            s[7]        u[1]       u[2]
#> s[1]  0.0000000 -0.89024499  3.9272674
#> s[2]  0.0000000  1.16133026  0.3454873
#> s[3]  0.0000000 -2.03050494 -0.3256890
#> s[4]  1.0000000 -0.93037412 -0.5852534
#> s[5]  1.6063890  0.01281744 -1.0897144
#> s[6] -0.2244914 -0.35879062  0.7272370
#> s[7] -2.2560126 -1.02087405  2.5244620
#> x[1]  0.0000000  1.00000000  0.0000000
#> x[2]  0.0000000  0.00000000  1.0000000
#> x[3]  0.0000000  0.00000000  0.0000000
#> 
#> $Hsv
#> NULL
#> 
#> $nu
#> [1] 3 2 2
#> 
# check impulse response
all.equal(pseries(out$Xs, lag.max = 20), Xi)
#> [1] TRUE

Xs1 = as.stsp(Xi, method = 'echelon')
all.equal(Xs1, out$Xs)
#> [1] TRUE
```
