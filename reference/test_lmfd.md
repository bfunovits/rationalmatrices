# Create Test Rational Matrix in LMFD Form

This simple tool may be used to create a random, \\(m,n)\\-dimensional,
rational matrix in LMFD form \$\$x(z)=a^{-1}(z) b(z)\$\$ The degrees of
the polynomials \\a(z), b(z)\\ is denoted with \\p\\ and \\q\\
respectively.

## Usage

``` r
test_lmfd(
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

[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
object, which represents the generated rational matrix \\x(z)\\ in LMFD
form.

## Details

We require \\m\>0\\ and \\p\geq 0\\. The left factor \\a(z)\\ is
normalized as \\a(0)=I_m\\ where \\I_m\\ denotes the
\\(m,m)\\-dimensional identity matrix.

The user may prescribe lower bounds for the moduli of the zeroes and/or
poles of the rational matrix. In this case the procedure simply
generates (up to n.trials) random matrices until a matrix is found which
satisfies the constraint. The standard deviation of the normal
distribution, which is used to generate the random entries, is decreased
in each step. Of course this is a very crude method and it may fail or
need a very large number of randomly generated matrices.

## Examples

``` r
### generate a random (2 x 2) rational matrix in LMFD form with degrees p=1 and q =1
### we require that the matrix has no poles and no zeroes within the unit circle!
x = try(test_lmfd(dim = c(2,2), degrees = c(1,1), digits = 2, bpoles = 1, bzeroes = 1))
if (!inherits(x, 'try-error')) {
   print(x)
   print(abs(poles(x)))
   print(abs(zeroes(x)))
}
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     0    -0.56 -0.53
#> [2,]        0     1     0.70 -1.08
#> right factor b(z):
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]     0.17  0.66     0.02 -0.75
#> [2,]     0.35 -0.25     0.27 -0.42
#> [1] 1.012324 1.012324
#> [1] 1.166866 1.207566
```
