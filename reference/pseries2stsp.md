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
#>              s[1]         s[2]        s[3]        s[4]        s[5]        s[6]
#> s[1]  0.183499666  0.728423835 -0.03469093  0.30605270 -0.24882405 -0.03689276
#> s[2] -0.884307257 -0.047354629 -0.24399663  0.22769285 -0.16287316 -0.02231447
#> s[3]  0.148429873 -0.566354944 -0.09814646  0.11074722 -0.14025660 -0.27280450
#> s[4] -0.149364712 -0.035793461  0.80868358  0.42419863  0.19278394 -0.08510797
#> s[5] -0.166325711  0.056117924  0.19029604 -0.62918729 -0.18612023  0.07306016
#> s[6] -0.036395496  0.065231455 -0.16983919  0.01824769  0.34901112 -0.67276767
#> s[7] -0.002712822  0.004381256 -0.04606908  0.05136175 -0.05207941  0.09274882
#> x[1]  0.695696769 -1.654240318  0.63773351  0.72097874 -1.45449842  0.36363193
#> x[2] -1.426338738  0.123081058  0.68972060  0.05854797  0.57146318  0.50458473
#> x[3]  0.853864196 -0.652161618 -1.21425314  1.42671842  0.96958912  0.58640875
#>             s[7]        u[1]        u[2]
#> s[1] -0.02847378  0.63654436  2.42256315
#> s[2]  0.01339733  0.98293014  0.48929830
#> s[3] -0.04920701 -0.03861519  1.95765533
#> s[4]  0.01682560 -0.75708394  0.19314649
#> s[5] -0.07975451 -0.94704506  0.72419834
#> s[6] -0.06531945 -0.68669872  0.10011843
#> s[7]  0.15114895 -0.28441271 -0.01664408
#> x[1] -0.04268630  1.00000000  0.00000000
#> x[2] -0.24789978  0.00000000  1.00000000
#> x[3] -0.10822461  0.00000000  0.00000000
#> 
#> $Hsv
#>  [1] 1.489987e+01 1.376184e+01 9.013358e+00 8.420514e+00 5.660010e+00
#>  [6] 2.735623e+00 1.635695e-01 1.915554e-15 1.403902e-15 1.165775e-15
#> [11] 8.833070e-16 8.329450e-16 5.852340e-16 5.483620e-16 3.897172e-16
#> [16] 3.643068e-16 3.043300e-16 2.425626e-16 1.915429e-16 1.703287e-16
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
#>            s[1]       s[2]       s[3]       s[4]       s[5]      s[6]
#> s[1]  0.0000000  0.0000000  0.0000000  1.0000000  0.0000000  0.000000
#> s[2]  0.0000000  0.0000000  0.0000000  0.0000000  1.0000000  0.000000
#> s[3]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  1.000000
#> s[4]  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> s[5]  0.7184332 -2.2054399  1.6221187 -2.1820028  0.9047551  3.213329
#> s[6] -0.5011707  0.1645196 -0.1596492  0.1252769  1.0032481 -1.469158
#> s[7] -0.3999640  2.1731533 -1.5464019  2.3107985 -0.7536741 -3.044312
#> x[1]  1.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> x[2]  0.0000000  1.0000000  0.0000000  0.0000000  0.0000000  0.000000
#> x[3]  0.0000000  0.0000000  1.0000000  0.0000000  0.0000000  0.000000
#>            s[7]       u[1]       u[2]
#> s[1]  0.0000000 -0.6137179  1.2474407
#> s[2]  0.0000000 -1.6751004 -1.5651318
#> s[3]  0.0000000 -2.4209141  0.4106014
#> s[4]  1.0000000  0.2821515  6.1027806
#> s[5] -0.6361931 -1.0575274 -1.4841053
#> s[6]  1.7541334  1.1874504  4.2300183
#> s[7]  0.3188609  0.8803605  2.1762808
#> x[1]  0.0000000  1.0000000  0.0000000
#> x[2]  0.0000000  0.0000000  1.0000000
#> x[3]  0.0000000  0.0000000  0.0000000
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
