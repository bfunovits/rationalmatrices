# Jacobian of the Solution of the Lyapunov Equation

This (internal helper) function considers the solution of a Lyapunov
equation \$\$P = A P A' + Q\$\$ where \\A,Q\\ are real valued, square
matrices and \\Q\\ is symmetric. The *directional derivative* of \\P\\
along (the matrices) \\dA\\, \\dQ\\ is given by the solution of the
Lyapunov Equation: \$\$dP = A dP A' + dA P A' + A P dA' + dQ\$\$

## Usage

``` r
lyapunov_Jacobian(A, Q, dA, dQ, non_stable = c("ignore", "warn", "stop"))
```

## Arguments

- A, Q:

  \\(m,m)\\ matrices. Note that the routine silently assumes that \\Q\\
  is symmetric (and hence the solution \\P\\ is also symmetric).

- dA, dQ:

  \\(m^2,n)\\ matrices. Each column of \\dA\\, \\dQ\\ determines a
  direction along which the derivative is computed. Note that the
  routine silently assumes that each column of \\dQ\\ represents a
  symmetric matrix (and hence \\dP\\ is also symmetric).

- non_stable:

  (character string) indicates what to do, when \\A\\ is not stable.

## Value

List with the slots

- P:

  (\\(m^2,n)\\-dimensional matrix) Solution of the lyapunov equation

- J:

  (\\(m^2,n)\\-dimensional matrix) Jacobian of the vectorised solution
  of the Lyapunov equation. Each column of \\J\\ is the directional
  derivative of \\vec(P)\\ along the respective columns of \\dA\\ and
  \\dQ\\.

- lambda:

  Eigenvalues of \\A\\.

- is_stable:

  (boolean) Is \\A\\ stable or not?

## See also

[`lyapunov`](https://bfunovits.github.io/rationalmatrices/reference/lyapunov.md).

## Examples

``` r
m = 5

A = matrix(rnorm(m^2), nrow = m, ncol = m)
Q = crossprod(matrix(rnorm(m^2), nrow = m, ncol = m))
P = lyapunov(A, Q)
all.equal(Q, P - A %*% P %*% t(A))
#> [1] TRUE

n = 6
dA = matrix(rnorm((m^2)*n), nrow = m^2, ncol = n)
dQ = matrix(rnorm((m^2)*n), nrow = m^2, ncol = n)
for (i in (1:n)) {
  # make sure that each column of dQ corresponds to a symmetric matrix!
  junk = matrix(dQ[,i], nrow = m, ncol = m)
  junk = junk + t(junk)
  dQ[,i] = junk
}
out = lyapunov_Jacobian(A, Q, dA, dQ)
all.equal(out$P, P)
#> [1] TRUE
 
eps = 1e-8
theta = rnorm(n)
matrix(out$J %*% theta, nrow = m, ncol = m)
#>            [,1]       [,2]      [,3]      [,4]      [,5]
#> [1,]   73.79928 -133.42202 -91.52057 -36.99219 -69.90643
#> [2,] -133.42202  248.21825 144.23692  63.61531 107.13161
#> [3,]  -91.52057  144.23692 347.25351 -23.82065 324.63611
#> [4,]  -36.99219   63.61531 -23.82065  57.38041 -33.49289
#> [5,]  -69.90643  107.13161 324.63611 -33.49289 290.65510

# compute the derivative via "finite differences"
dP = lyapunov(A + matrix(dA %*% theta, nrow = m, ncol = m)*eps, 
              Q + matrix(dQ %*% theta, nrow = m, ncol = m)*eps)
all.equal(matrix(out$J %*% theta, nrow = m, ncol = m), 
          (dP - P)/eps, scale = mean(abs(out$J)), tol = 1e-6)
#> [1] TRUE
```
