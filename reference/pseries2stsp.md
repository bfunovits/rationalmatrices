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
#> s[1] -0.918255846  0.11314696 -0.11434858 -0.01966634 -0.03426064 -0.04574571
#> s[2] -0.262942985  0.22248080  0.17660486  0.48595070 -0.09831646  0.07716176
#> s[3] -0.145569035 -0.56562997 -0.04978727 -0.35366320 -0.29539687  0.14677140
#> s[4] -0.018960484  0.10368292  0.88873284 -0.08031917 -0.06741444 -0.08717843
#> s[5] -0.022424424  0.19341776 -0.11719190 -0.35868311  0.47789770 -0.15488296
#> s[6] -0.006242367 -0.08155169  0.04585152  0.09001492  0.44649802  0.35380768
#> s[7]  0.007638275 -0.02923697  0.09236008 -0.14623327  0.01393301 -0.12807702
#> x[1] -1.223252078 -0.07946179  0.76095605 -1.02079802  0.36643964  0.25826828
#> x[2] -0.679048680 -1.15913727  0.16323898  0.77932872  0.66638605  0.21023515
#> x[3] -0.340243981 -2.20486302  0.18979491  0.39746228 -0.06532628 -0.31718086
#>             s[7]        u[1]       u[2]
#> s[1]  0.02690261 -0.72302169 -1.4621133
#> s[2] -0.01435545  1.90957858  1.4686526
#> s[3]  0.11006682 -0.09821544  1.0258142
#> s[4]  0.09538677 -0.77132241  0.3937009
#> s[5]  0.20238877 -0.07683888  0.4591242
#> s[6]  0.35323813  0.02733752 -0.0684763
#> s[7]  0.02836679  0.39641318 -0.2876255
#> x[1] -0.35193988  1.00000000  0.0000000
#> x[2]  0.10223858  0.00000000  1.0000000
#> x[3] -0.12782394  0.00000000  0.0000000
#> 
#> $Hsv
#>  [1] 1.496562e+01 8.475439e+00 4.820186e+00 4.662340e+00 1.597350e+00
#>  [6] 5.450915e-01 3.960418e-01 8.422564e-16 7.246733e-16 6.396904e-16
#> [11] 4.525230e-16 3.823088e-16 3.325051e-16 2.929219e-16 2.519187e-16
#> [16] 1.921928e-16 1.625905e-16 1.399919e-16 1.268113e-16 1.109006e-16
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
#>            s[1]      s[2]      s[3]      s[4]       s[5]      s[6]       s[7]
#> s[1]  0.0000000  0.000000  0.000000  1.000000  0.0000000  0.000000  0.0000000
#> s[2]  0.0000000  0.000000  0.000000  0.000000  1.0000000  0.000000  0.0000000
#> s[3]  0.0000000  0.000000  0.000000  0.000000  0.0000000  1.000000  0.0000000
#> s[4]  0.0000000  0.000000  0.000000  0.000000  0.0000000  0.000000  1.0000000
#> s[5]  1.3657737 -4.593689  3.287473  2.674436  2.3652899 -5.338705  2.1100132
#> s[6]  0.5070613 -1.784064  1.274413  1.253389  0.1357087 -1.486196  1.3810186
#> s[7] -1.3143465  3.942287 -2.825066 -2.147438 -2.3066104  5.046795 -0.8449032
#> x[1]  1.0000000  0.000000  0.000000  0.000000  0.0000000  0.000000  0.0000000
#> x[2]  0.0000000  1.000000  0.000000  0.000000  0.0000000  0.000000  0.0000000
#> x[3]  0.0000000  0.000000  1.000000  0.000000  0.0000000  0.000000  0.0000000
#>            u[1]        u[2]
#> s[1]  1.2847160  2.30232470
#> s[2] -2.3445714  0.02690289
#> s[3] -4.3438902 -2.36103502
#> s[4] -1.6345174 -3.42879544
#> s[5] -0.3548682 -1.34275294
#> s[6] -0.8592594 -2.56748588
#> s[7]  1.4037512  1.42852462
#> x[1]  1.0000000  0.00000000
#> x[2]  0.0000000  1.00000000
#> x[3]  0.0000000  0.00000000
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
