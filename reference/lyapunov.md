# Lyapunov Equation

This function solves the Lyapunov equation \$\$P = A P A' + Q\$\$ where
\\A,Q\\ are real valued, square matrices and \\Q\\ is symmetric. The
Lyapunov equation has a unique solution if \\\lambda_i(A)\lambda_j(A)
\neq 1\\ holds for all eigenvalues of \\A\\. If \\A\\ is stable (i.e.
the spectral radius of \\A\\ is less than one) and \\Q\\ is positive
semidefinite then the solution \\P\\ is also positive semidefinite.  
The procedure uses the Schur decomposition(s) of \\A\\ and computes the
solution "column by column", see (Kitagawa 1977; Hammarling 1982) .

## Usage

``` r
lyapunov(A, Q, non_stable = c("ignore", "warn", "stop"), attach_lambda = FALSE)
```

## Arguments

- A, Q:

  \\(m,m)\\ matrices. Note that the routine silently assumes that \\Q\\
  is symmetric (and hence the solution \\P\\ is also symmetric).

- non_stable:

  (character string) indicates what to do, when \\A\\ is not stable.

- attach_lambda:

  (boolean) if yes, then the eigenvalues of \\A\\ are attached to the
  solution \\P\\ as an attribute.

## Value

P (\\(m,m)\\ matrix) the solution of the Lyapunov equation.

## References

Kitagawa G (1977). “An algorithm for solving the matrix equation X = FXF
T + S.” *International Journal of Control*, **25**(5), 745-753.
[doi:10.1080/00207177708922266](https://doi.org/10.1080/00207177708922266)
, http://dx.doi.org/10.1080/00207177708922266,
<http://dx.doi.org/10.1080/00207177708922266>.

Hammarling SJ (1982). “Numerical solution of the stable, nonnegative
definite Lyapunov equation.” *IMA Journal of Numerical Analysis*,
**2**(3), 303–323. ISSN 0272-4979 (print), 1464-3642 (electronic).

## Examples

``` r
# A is stable and Q is positve definite
m = 4
A = diag(runif(m, min = -0.9, max = 0.9))
V = matrix(rnorm(m*m), nrow = m, ncol = m)
A = V %*% A %*% solve(V)
B = matrix(rnorm(m*m), nrow = m, ncol = m)
Q = B %*% t(B)
P = lyapunov(A, Q)
all.equal(P, A %*% P %*% t(A) + Q)
#> [1] TRUE

# unstable matrix A
A = diag(runif(m, min = -0.9, max = 0.9))
A[1,1] = 2
V = matrix(rnorm(m*m), nrow = m, ncol = m)
A = V %*% A %*% solve(V)
P = lyapunov(A, Q)
all.equal(P, A %*% P %*% t(A) + Q)
#> [1] TRUE
# note that the solution P (in general) is not positive semidefinite
eigen(P, only.values= TRUE, symmetric = TRUE)$values
#> [1]  23.0279705   4.2563159   0.7768525 -22.4080598

# attach the eigenvalues of A to the solution P
P = lyapunov(A, Q, attach = TRUE)
print(P)
#>           [,1]      [,2]       [,3]       [,4]
#> [1,]  2.252003  4.034818 -1.5223909 -2.3998956
#> [2,]  4.034818  3.675940 18.8459997 -9.8545848
#> [3,] -1.522391 18.846000 -6.3919721 -0.7017767
#> [4,] -2.399896 -9.854585 -0.7017767  6.1171076
#> attr(,"lambda")
#> [1]  2.0000000+0i -0.1898799+0i -0.7036820+0i -0.4964202+0i

# issue a warning message
P = lyapunov(A, Q, non_stable = 'warn')
#> Warning: "A" matrix is not stable

if (FALSE) { # \dontrun{
# throw an error
P = lyapunov(A, Q, non_stable = 'stop')
} # }
```
