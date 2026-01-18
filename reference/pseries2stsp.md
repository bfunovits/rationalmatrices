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
#>              s[1]        s[2]        s[3]        s[4]        s[5]        s[6]
#> s[1] -0.807300474 -0.32219180 -0.06313465 -0.11938096  0.08242245 -0.01442464
#> s[2] -0.434140586  0.22754640  0.19243253 -0.23548484 -0.15501905  0.08960222
#> s[3] -0.043127924 -0.72232909 -0.29146496  0.24321060 -0.22461734 -0.00050370
#> s[4]  0.130307806 -0.20736416  0.03781880 -0.25313031 -0.13904808  0.00143141
#> s[5]  0.002708488 -0.06158707  0.32110301  0.37461404 -0.16240254  0.15790083
#> s[6] -0.001019213  0.02806528 -0.04614370  0.00429105 -0.44472816  0.47266134
#> s[7] -0.003515291  0.01333606 -0.03077989  0.07008102 -0.04480428 -0.23256821
#> x[1] -0.547579641 -0.23092261  1.00074638  0.78111035  0.31877616  0.16804938
#> x[2] -1.372711562  1.26868467 -0.02259687  0.57779778 -0.32317167 -0.18959585
#> x[3]  0.278928295 -0.65438460  1.66153227 -0.20271246 -0.21398345 -0.18953977
#>              s[7]        u[1]        u[2]
#> s[1]  0.005924672  0.25515410  1.64821756
#> s[2] -0.006333880 -1.18142553 -1.43275009
#> s[3] -0.004538908  0.21595785 -0.97014138
#> s[4] -0.098935441 -1.07485543  0.46272738
#> s[5]  0.073161007 -0.01470511  0.04452475
#> s[6]  0.082539367  0.01581244  0.03123066
#> s[7] -0.239933093 -0.04216362  0.01548928
#> x[1] -0.063458265  1.00000000  0.00000000
#> x[2] -0.011404151  0.00000000  1.00000000
#> x[3]  0.034300644  0.00000000  0.00000000
#> 
#> $Hsv
#>  [1] 9.039049e+00 5.720047e+00 4.473838e+00 1.906312e+00 7.768727e-01
#>  [6] 2.176221e-01 3.175858e-02 5.080679e-16 3.625040e-16 2.636120e-16
#> [11] 2.281664e-16 1.972415e-16 1.565247e-16 1.089321e-16 8.850786e-17
#> [16] 6.828421e-17 6.536472e-17 5.917662e-17 5.186065e-17 3.800573e-17
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
#>             s[1]       s[2]       s[3]      s[4]       s[5]       s[6]
#> s[1]  0.00000000  0.0000000  0.0000000  1.000000  0.0000000  0.0000000
#> s[2]  0.00000000  0.0000000  0.0000000  0.000000  1.0000000  0.0000000
#> s[3]  0.00000000  0.0000000  0.0000000  0.000000  0.0000000  1.0000000
#> s[4]  0.00000000  0.0000000  0.0000000  0.000000  0.0000000  0.0000000
#> s[5] -0.20688490  0.6327729  0.3361774  1.509615 -1.3372046 -1.4775933
#> s[6]  0.23997806 -0.8847923 -0.3946407 -2.638055  1.7052262  2.1778571
#> s[7]  0.05088021 -0.3537938 -0.1616041 -1.053335  0.3826484  0.7365648
#> x[1]  1.00000000  0.0000000  0.0000000  0.000000  0.0000000  0.0000000
#> x[2]  0.00000000  1.0000000  0.0000000  0.000000  0.0000000  0.0000000
#> x[3]  0.00000000  0.0000000  1.0000000  0.000000  0.0000000  0.0000000
#>           s[7]        u[1]       u[2]
#> s[1]  0.000000 -0.48971573 -1.1626415
#> s[2]  0.000000 -2.47280335 -3.8114373
#> s[3]  0.000000  1.41968722 -0.3233334
#> s[4]  1.000000  0.74143118  2.3880109
#> s[5]  2.967802 -0.09153545 -0.3352290
#> s[6] -3.321830  0.94441383  2.8190945
#> s[7] -1.894676  0.42152157  0.4074680
#> x[1]  0.000000  1.00000000  0.0000000
#> x[2]  0.000000  0.00000000  1.0000000
#> x[3]  0.000000  0.00000000  0.0000000
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
