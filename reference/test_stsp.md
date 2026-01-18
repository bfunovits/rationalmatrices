# Create Test Rational Matrix in Statespace Form

This simple tool may be used to create a random, \\(m,n)\\-dimensional,
rational matrix in statespace form with statespace dimension \\s\\.

## Usage

``` r
test_stsp(
  dim = c(1, 1),
  s = NULL,
  nu = NULL,
  D = NULL,
  digits = NULL,
  bpoles = NULL,
  bzeroes = NULL,
  n.trials = 100
)
```

## Arguments

- dim:

  integer vector `c(m,n)`.

- s:

  state dimension (or NULL).

- nu:

  vector with the Kronecker indices (or `NULL`). Either the statespace
  dimension `s` or the Kronecker indices `nu` must be non `NULL`. If
  both parameters are given, then the parameter `s` is ignored.

- D:

  \\(m,n)\\ dimensional matrix (or `NULL`). See the details below.

- digits:

  integer, if non NULL then the randomly generated numbers are rounded
  to "digits" number of decimal places.

- bpoles:

  lower bound for the moduli of the poles of the rational matrix (or
  NULL).

- bzeroes:

  lower bound for the moduli of the zeroes of the rational matrix (or
  NULL). This parameter is ignored for non-square matrices (m != n).

- n.trials:

  maximum number of trials.

## Value

[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
object, which represents the generated rational matrix \\x(z)\\ in
statespace form.

## Details

If the Kronecker indices (parameter `nu`) are given, then a statespace
model in echelon canonical form is generated. This means that some of
the entries of the \\A,B,C\\ matrices are fixed to be one or zero and
the others are considerd as "free". See also
[`Kronecker-Indices`](https://bfunovits.github.io/rationalmatrices/reference/Kronecker-Indices.md).
The entries of the \\A, B, C\\ matrices, which are not a priori fixed
are randomly generated.

If only the state dimension \\s\\ (parameter `s`) is given, then all
entries of the \\A, B, C\\ matrices are considered as "free".

The \\D\\ matrix defaults to a \\(m,n)\\-dimensional diagonal matrix
with ones on the diagonal (`diag(x=1, nrow = m, ncol = n)`). However,
one may also pass an arbitray (compatible) \\D\\ matrix to the
procedure. This matrix may contain `NA`'s, which then are replaced by
random numbers.

The user may prescribe lower bounds for the moduli of the poles and the
zeroes of the rational matrix. In this case the procedure simply
generates (up to `n.trials`) random matrices until a matrix is found
which satisfies the constraints. The standard deviation of the normal
distribution, which is used to generate the random entries, is decreased
in each step. Of course this is a very crude method and it may fail or
need a very large number of randomly generated matrices.

Note also, that the generated model may be non-minimal.

## Examples

``` r
## random (2 x 2) statespace realization with state dimension s = 3
## no poles and zeroes within the unit circle
x = test_stsp(dim = c(2,2), s = 3, digits = 2, bpoles = 1, bzeroes = 1)
x
#> statespace realization [2,2] with s = 3 states
#>       s[1]  s[2]  s[3]  u[1]  u[2]
#> s[1] -0.22 -1.00 -0.02  1.06 -0.84
#> s[2]  0.53 -0.52 -0.21  0.24  0.33
#> s[3] -0.09  0.41  0.13 -0.40  0.86
#> x[1] -0.07 -0.48  0.09  1.00  0.00
#> x[2]  0.41  0.14 -0.13  0.00  1.00
min(abs(poles(x))) > 1
#> [1] TRUE
min(abs(zeroes(x))) > 1
#> [1] TRUE
pseries2nu(pseries(x, lag.max = 5)) # Kronecker indices
#> [1] 2 1

## random (3 x 2) statespace realization in echelon canonical form
## D is lower triangular (with ones on the diagonal)
## no poles within the unit circle
x = test_stsp(dim = c(3, 2), nu = c(2,3,0), D = matrix(c(1,NA,NA,0,1,NA), nrow = 3, ncol = 2), 
              digits = 2, bpoles = 1)

x
#> statespace realization [3,2] with s = 5 states
#>       s[1]  s[2]  s[3] s[4] s[5]  u[1]  u[2]
#> s[1]  0.00  0.00  1.00 0.00 0.00 -0.20 -1.51
#> s[2]  0.00  0.00  0.00 1.00 0.00  0.45  0.25
#> s[3] -0.39  0.24  0.00 0.03 0.00  1.00  0.23
#> s[4]  0.00  0.00  0.00 0.00 1.00  0.02  0.47
#> s[5] -0.31 -0.66 -0.23 0.05 0.22  0.97 -0.92
#> x[1]  1.00  0.00  0.00 0.00 0.00  1.00  0.00
#> x[2]  0.00  1.00  0.00 0.00 0.00  0.01  1.00
#> x[3]  0.60  0.06  0.00 0.00 0.00 -0.26  1.17
min(abs(poles(x))) > 1
#> [1] TRUE
pseries2nu(pseries(x, lag.max = 10)) # check Kronecker indices
#> [1] 2 3 0
```
