# Create Test Rational Matrix in RMFD Form

This simple tool may be used to create a random, \\(m,n)\\-dimensional,
rational matrix in RMFD form \$\$k(z) = d(z) c^{-1}(z)\$\$ The degrees
of the polynomials \\(c'(z), d'(z))'\\ is denoted with \\p\\ and \\q\\
respectively.

## Usage

``` r
test_rmfd(
  dim = c(1, 1),
  degrees = c(1, 1),
  digits = NULL,
  bpoles = NULL,
  bzeroes = NULL,
  n.trials = 100
)
```

## Arguments

- dim:

  integer vector `c(m,n)`.

- degrees:

  integer vector `c(p,q)`.

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

[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
object, which represents the generated rational matrix \\k(z)\\ in RMFD
form.

## Details

We require \\m\>0\\ and \\p\geq 0\\. The right factor \\c(z)\\ is
normalized as \\c(0)=I_n\\ where \\I_n\\ denotes the
\\(n,n)\\-dimensional identity matrix.

The user may prescribe lower bounds for the moduli of the zeroes and/or
poles of the rational matrix. In this case, the procedure simply
generates (up to n.trials) random matrices until a matrix is found which
satisfies the constraint. The standard deviation of the normal
distribution, which is used to generate the random entries, is decreased
in each step. Of course, this is a very crude method and it may fail or
need a very large number of randomly generated matrices.

## Examples

``` r
### generate a random (2 x 2) rational matrix in RMFD form with degrees p=1 and q=1
### we require that the matrix has no poles and no zeroes within the unit circle!
x = try(test_rmfd(dim = c(2,2), degrees = c(1,1), digits = 2, bpoles = 1, bzeroes = 1))
if (!inherits(x, 'try-error')) {
   print(x)
   print(abs(poles(x)))
   print(abs(zeroes(x)))
}
#> ( 2 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]     0.32  0.23    -0.28  0.43
#> [2,]    -0.49  0.49     0.44 -0.22
#> right factor c(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0     0.82  0.36
#> [2,]        0     1    -0.07 -0.04
#> [1]   1.266424 103.898003
#> [1] 1.118871 1.887679
```
