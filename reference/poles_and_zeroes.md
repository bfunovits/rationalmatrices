# Poles and Zeroes

Compute the poles and zeroes of a rational matrix. For polynomial
matrices and rational matrices in left matrix fraction form the poles
(and zeroes) are computed via the (reciprocals of the) eigenvalues of
the associated companion matrices, see also
[`companion_matrix`](https://bfunovits.github.io/rationalmatrices/reference/companion_matrix.md).
For statespace realizations the poles are computed via (the reciprocals
of) the eigenvalues of the state transition matrix \\A\\ and for the
zeroes the eigenvalues of the state transition matrix \\A-BD^{-1}C\\ of
the inverse are used.

## Usage

``` r
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'polm'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'lpolm'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'lmfd'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'rmfd'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'stsp'
zeroes(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'polm'
poles(x, ...)

# S3 method for class 'lmfd'
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'rmfd'
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)

# S3 method for class 'stsp'
poles(x, tol = sqrt(.Machine$double.eps), print_message = TRUE, ...)
```

## Arguments

- x:

  an object which represents a rational matrix (i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  or
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object).

- tol:

  Double. Default set to `sqrt(.Machine$double.eps)`. Required to decide
  on when a root is to be considered "at infinity".

- print_message:

  Boolean. Default set to TRUE. Prints a message if roots "at infinity "
  are discarded.

- ...:

  not used.

## Value

Vector of poles, respectively zeroes.

## Details

The methods do not return numerically reliable and correct results in
all cases. For some more details see the vignette [Rational
Matrices](https://bfunovits.github.io/rationalmatrices/doc/rational_matrices.md).

- Zeroes are only computed for square, non singular matrices which have
  no zero at \\z=0\\. If the matrix evaluated at \\z=0\\ is close to
  singular, the results may be unreliable.

- The procedures use a threshold `tol` in order to decide whether a
  small eigenvalue returned by
  [`eigen`](https://rdrr.io/r/base/eigen.html) corresponds to a "true
  zero" eigenvalue or not.

- If the pair \\a,b\\ of polynomials of the LMFD is not left coprime
  then a pole/zero cancellation occurs. This is not taken into account
  by the procedures. Hence, in this case, the results also contain some
  spurious poles/zeroes. This happens also for non minimal state space
  realizations.

## Examples

``` r
 
# zeroes of polynomial matrices #############################################

# scalar polynomial ###
(a = polm(c(1, 0, 0, 0.5, 0)))
#> ( 1 x 1 ) matrix polynomial with degree <= 4 
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1] z^4 [,1]
#> [1,]        1        0        0      0.5        0
(z = zeroes(a))
#> [1] -1.2599210+0.000000i  0.6299605-1.091124i  0.6299605+1.091124i

# compare with the result of "polyroot"
all.equal(sort(z), sort(polyroot(as.vector(a))))
#> [1] TRUE

# zero degree polynomial (have no zeroes) ###
zeroes(polm(diag(3)))
#> numeric(0)

# (2 x 2) polynomial of degree 2 ### 
a = polm(dbind(d = 3, diag(2), test_array(dim = c(2,2,2))))
(z = zeroes(a))
#> [1] -0.002994225  0.265114779 -0.991220211 -1.270900343

# check the rank of a(z) at the computed zeroes 
az = zvalues(a, z)
apply(az, MARGIN = 3, FUN = function(x) {d = svd(x)$d; min(d)/max(d)})
#> [1] 1.403359e-16 8.297615e-17 1.374392e-14 0.000000e+00

if (FALSE) { # \dontrun{
# the following examples throw an error
zeroes(polm(c(0, 0, 0, 0.5))) # constant term is zero
zeroes(polm(test_array(dim = c(2, 1, 3)))) # non-square polynomial
zeroes(polm(test_array(dim = c(2, 2, 0)))) # zero polynomial
} # }
 
# zeroes of a Laurent polynomial #################################
(lp = lpolm(1:5, min_deg = -7))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= -3, and minimal degree >= -7
#>      z^-7 [,1] z^-6 [,1] z^-5 [,1] z^-4 [,1] z^-3 [,1]
#> [1,]         1         2         3         4         5
(p = polm(1:5))
#> ( 1 x 1 ) matrix polynomial with degree <= 4 
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1] z^4 [,1]
#> [1,]        1        2        3        4        5
zeroes(p)
#> [1] -0.5378323-0.3582847i -0.5378323+0.3582847i  0.1378323-0.6781544i
#> [4]  0.1378323+0.6781544i
zeroes(lp)
#> [1] -0.5378323-0.3582847i -0.5378323+0.3582847i  0.1378323-0.6781544i
#> [4]  0.1378323+0.6781544i
 
# zeroes of a rational matrix in LMFD form #################################

c = lmfd(test_polm(dim = c(2,2), degree = 3, random = TRUE),
         test_polm(dim = c(2,2), degree = 1, random = TRUE))
(z = zeroes(c))
#> [1] -0.3875626-0.4351629i -0.3875626+0.4351629i
all.equal(z, zeroes(c$b))
#> [1] TRUE
 
# zeroes of a rational matrix in RMFD form #################################

k = rmfd(c = test_polm(dim = c(2,2), degree = 3, random = TRUE),
         d = test_polm(dim = c(2,2), degree = 1, random = TRUE))
(z = zeroes(k))
#> [1]  1.654498 -2.717882
all.equal(z, zeroes(k$d))
#> [1] TRUE
 
# zeroes of a rational matrix in statespace form ###########################

k = stsp(A = matrix(rnorm(3*3), nrow = 3, ncol = 3),
         B = matrix(rnorm(3*2), nrow = 3, ncol = 2),
         C = matrix(rnorm(3*2), nrow = 2, ncol = 3),
         D = matrix(rnorm(2*2), nrow = 2, ncol = 2))
(z = zeroes(k, tol = 0))
#> [1] -0.05182903-0.3429019i -0.05182903+0.3429019i 24.69212861+0.0000000i
all.equal(z, 1/(eigen(k$A - k$B %*% solve(k$D, k$C), only.values = TRUE)$values))
#> [1] TRUE

if (FALSE) { # \dontrun{
k = stsp(k$A, k$B, k$C, 
         D = matrix(rnorm(2*1), nrow = 2, ncol = 1)[,c(1,1)])  # D is singular
zeroes(k)                                                      # zeroes() throws an error

k = stsp(k$A, k$B[,1,drop = FALSE], k$C, k$D[,1,drop = FALSE]) # (2 x 1) rational matrix
zeroes(k)                          # throws an error, since k is not square
} # }

# poles of polynomial matrices #############################################

# polynomials have no poles ###
poles(test_polm(dim = c(2,1), degree = 2, random = TRUE)) 
#> numeric(0)

# poles of a rational matrix in LMFD form ##################################

(z = poles(c))
#> [1]  0.1227278+0.000000i  0.1960404+0.000000i  0.5171242-1.309299i
#> [4]  0.5171242+1.309299i -1.0682088-1.499647i -1.0682088+1.499647i
all.equal(z, zeroes(c$a))
#> [1] TRUE

# poles of a rational matrix in RMFD form ##################################

(z = poles(c))
#> [1]  0.1227278+0.000000i  0.1960404+0.000000i  0.5171242-1.309299i
#> [4]  0.5171242+1.309299i -1.0682088-1.499647i -1.0682088+1.499647i
all.equal(z, zeroes(c$a))
#> [1] TRUE

# poles of a rational matrix in statespace form ###########################

(z = poles(k, tol = 0))
#> [1] -0.5266464-0.4236672i -0.5266464+0.4236672i  1.1797597+0.0000000i
all.equal(z, 1/(eigen(k$A, only.values = TRUE)$values))
#> [1] TRUE
```
