# Hermitean Transpose

The *Hermitean transpose* of a rational matrix \\x(z)\\ is defined as
\$\$x^{\*}(z)=\overline{x(\bar{z}^{-1})}'.\$\$ This means e.g. for a
polynomial with real coefficients \\x(z)=a_0 + a_1 z +\cdots + a_p z^p\\
that the coefficient matrices are transposed and that \\z\\ is replaced
by \\z^{-1}\\: \\x^\*(z)=a_0' + a_1' z^{-1} +\cdots + a_p' z^{-p}\\,
i.e. the result is (in general) a Laurent polynomial.

## Usage

``` r
Ht(x)

# S3 method for class 'polm'
Ht(x)

# S3 method for class 'lpolm'
Ht(x)

# S3 method for class 'stsp'
Ht(x)

# S3 method for class 'zvalues'
Ht(x)
```

## Arguments

- x:

  rational matrix object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md).

## Value

A rational matrix object which represents the Hermitean transpose
\\x^\*(z)\\. The Hermitean transpose of a polynomial matrix (in general)
is a Laurent polynomial, therefore the output is of class `lpolm` if `x`
is a `polm` object. In all other cases the output is of the same class
as the input `x`.

## Details

The Hermitean transpose is only implemented for polynomial matrices
([`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
objects), Laurent polynomials
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
objects), rational matrices in statespace form
([`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
objects) and for the frequency response of rational matrices
([`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
objects).

The case of rational matrices in LMFD or RMFD form has not (yet) been
implemented. It is not possible to directly construct the Hermitean
transpose from given power series coefficients. Therefore `Ht()` does
not support
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
objects.

Finally note that the Hermitean transpose of a rational matrix
\\x(z)=D+zB(I-Az)^{-1}B\\ in statespace form (in general) has a pole at
zero, if the state transition matrix \\A\\ is singular. Since rational
matrices with a pole at \\z=0\\ have **no** statespace realization, the
procedure `Ht()` throws an error in this case.

## Examples

``` r
# rational matrix in statespace form 
x = test_stsp(dim = c(4,4), s = 2)
Ht(x)
#> statespace realization [4,4] with s = 2 states
#>            s[1]        s[2]         u[1]        u[2]       u[3]        u[4]
#> s[1]  0.3671372 -3.66425753  1.033860142  2.25709229 -1.8831275  -8.2356401
#> s[2] -0.3718365 -0.38385156  0.110588535 -0.01857762 -0.6684571  -0.1152447
#> x[1]  0.7893360 -5.77966931  2.629546389  3.69081328 -2.7288256 -13.3731873
#> x[2]  0.4898776 -7.02692920  1.983821925  5.19529112 -3.8572312 -15.4032514
#> x[3]  0.2212180  0.01568305 -0.005666152  0.12881495  1.2639137  -0.3706329
#> x[4]  0.5096678 -3.16970376  0.893253019  2.07184844 -1.4083783  -6.4740390
all.equal(zvalues(Ht(x)), Ht(zvalues(x)))
#> [1] TRUE

# Note (x(z)  x^*(z)) is Hermitean for |z| = 1
xx = zvalues(x) %r% Ht(zvalues(x))
# print(xx, digits = 2)
apply(unclass(xx), MARGIN = 3, FUN = isSymmetric)
#> [1] TRUE TRUE TRUE TRUE TRUE

# polynomial matrix 
x = test_polm(dim = c(2,3), degree = 1)
Ht(x)
#> ( 3 x 2 ) Laurent polynomial matrix with degree <= 0, and minimal degree >= -1
#>      z^-1 [,1]  [,2] z^0 [,1]  [,2]
#> [1,]       111   211      110   210
#> [2,]       121   221      120   220
#> [3,]       131   231      130   230

if (FALSE) { # \dontrun{
Ht(test_lmfd(dim = c(2,2), degrees = c(3,3)))
Ht(test_rmfd(dim = c(2,2), degrees = c(3,3)))
Ht(pseries(test_lmfd(dim = c(2,2), degrees = c(3,3))))
} # }
```
