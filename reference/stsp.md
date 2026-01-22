# Constructor for Statespace Realizations

Any rational (\\(m,n)\\-dimensional) matrix \\x(z)\\ which has no pole
at \\z=0\\ may be represented as \$\$x(z) = C(z^{-1} I_s - A)^{-1} B +
D\$\$ Here \\I_s\\ denotes the \\(s,s)\\-dimensional identity matrix and
\\A,B,C,D\\ are (real or complex valued) matrices of size \\(s,s)\\,
\\(s,n)\\, \\(m,s)\\ and \\(m,n)\\ respectively. The integer \\s \geq
0\\ is called the *state dimension* of the above *statespace
realization* of \\x(z)\\.

## Usage

``` r
stsp(A, B, C, D)
```

## Arguments

- A:

  \\(s,s)\\ matrix (or a vector of length \\(s^2)\\)

- B:

  \\(s,n)\\ matrix (or a vector of length \\(sn)\\)

- C:

  \\(m,s)\\ matrix (or a vector of length \\(ms)\\)

- D:

  \\(m,n)\\ matrix (or a vector of length \\(mn)\\)

## Value

An object of class `stsp`.

## Details

Internally statespace realizations are stored as an \\(s+m,
s+n)\\-dimensional matrix with an attribute `order = c(m,n,s)` and a
class attribute `c("stsp","ratm")`.

Any of the integers \\m,n,s\\ may be zero.

If the arguments `A,B,C` are missing, then a statespace realization with
statespace dimension \\s=0\\ is constructed. In this case `D` must a be
matrix or a scalar. If `A,B,C` are given (and compatible) and `D` is
missing, then `D = diag(x = 1, nrow = m, ncol = n)` is used.

## See also

- [`test_stsp`](https://bfunovits.github.io/rationalmatrices/reference/test_stsp.md)
  generates random polynomials.

- checks:
  [`is.stsp`](https://bfunovits.github.io/rationalmatrices/reference/is.md),
  [`is.miniphase`](https://bfunovits.github.io/rationalmatrices/reference/check.md),
  [`is.stable`](https://bfunovits.github.io/rationalmatrices/reference/check.md)
  and
  [`is.minimal`](https://bfunovits.github.io/rationalmatrices/reference/is.minimal.md).

- generic S3 methods: [`dim`](https://rdrr.io/r/base/dim.html),
  [`str`](https://rdrr.io/r/utils/str.html) and
  [`print`](https://rdrr.io/r/base/print.html).

- arithmetics:
  [`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md),
  matrix multiplication
  [`%r%`](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md),
  ...

- matrix operations: transpose
  [`t.stsp`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md),
  Hermitean transpose
  [`Ht`](https://bfunovits.github.io/rationalmatrices/reference/Ht.md),
  [`bind`](https://bfunovits.github.io/rationalmatrices/reference/bind.md),
  [`[.stsp`](https://bfunovits.github.io/rationalmatrices/reference/extract.md),
  ...

- extract the system matrices \\A,B,C,D\\ with
  [`$.stsp`](https://bfunovits.github.io/rationalmatrices/reference/extract.md).

- statespace realizations related tools: state transfromation
  [`state_trafo`](https://bfunovits.github.io/rationalmatrices/reference/state_trafo.md),
  observability and controllability matrices
  [`obs_matrix`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md)
  and
  [`ctr_matrix`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md),
  Grammian matrices
  [`grammians`](https://bfunovits.github.io/rationalmatrices/reference/grammians.md)
  and balanced realizations
  [`balance`](https://bfunovits.github.io/rationalmatrices/reference/balance.md)..

- [`reflect_poles`](https://bfunovits.github.io/rationalmatrices/reference/reflect_poles.md)
  and
  [`reflect_zeroes`](https://bfunovits.github.io/rationalmatrices/reference/reflect_zeroes.md)
  may be used to reflect poles and zeroes of a rational matrix in
  statespace form by multiplication with allpass rational matrices.

- [`zeroes`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md),
  ...

## Examples

``` r
### x(z) =  I
stsp(D = diag(3))
#> statespace realization [3,3] with s = 0 states
#>      u[1] u[2] u[3]
#> x[1]    1    0    0
#> x[2]    0    1    0
#> x[3]    0    0    1

### random (2 x 3) rational matrix in statespace form with state dimension s = 4
x = stsp(A = stats::rnorm(4*4), B = stats::rnorm(4*3), C = stats::rnorm(2*4))

is.stsp(x)
#> [1] TRUE
dim(x)
#> m n s 
#> 2 3 4 
str(x)
#> ( 2 x 3 ) statespace realization with s = 4 states
print(x, digits = 3)
#> statespace realization [2,3] with s = 4 states
#>        s[1]  s[2]   s[3]   s[4]   u[1]   u[2]   u[3]
#> s[1]  1.226 1.319 -0.139 -0.786  1.444 -0.178  1.714
#> s[2] -0.294 1.218  0.887 -1.421 -0.308 -1.376 -0.873
#> s[3] -1.177 0.794 -0.298 -0.403  0.312 -0.586 -1.099
#> s[4] -1.750 1.903 -0.739  0.901  0.031 -0.870  0.191
#> x[1] -0.055 1.591  0.140  0.415  1.000  0.000  0.000
#> x[2]  0.480 2.093  1.438 -1.124  0.000  1.000  0.000
```
