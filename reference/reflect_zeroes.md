# Reflect Zeroes of Rational Matrices

Given a square, rational matrix \\x(z)\\ and a set of zeroes of
\\x(z)\\, the function `reflect_zeroes` constructs a rational, allpass
matrix \\U(z)\\ such that \\x(z)U(z)\\ has zeroes which are the
reciprocals of the selected zeroes. The other zeroes of \\x(z)\\ and the
poles of \\x(z)\\ are not changed.

## Usage

``` r
reflect_zeroes(x, zeroes, ...)

# S3 method for class 'polm'
reflect_zeroes(
  x,
  zeroes,
  tol = sqrt(.Machine$double.eps),
  check_zeroes = TRUE,
  ...
)

# S3 method for class 'lmfd'
reflect_zeroes(
  x,
  zeroes,
  tol = sqrt(.Machine$double.eps),
  check_zeroes = TRUE,
  ...
)

# S3 method for class 'stsp'
reflect_zeroes(x, zeroes, tol = sqrt(.Machine$double.eps), ...)
```

## Arguments

- x:

  Square, rational matrix object (with real coefficients).

- zeroes:

  Complex or real vector of zeroes to be reflected. It is assumed that
  for complex conjugated pairs of roots, only **one** is contained in
  this vector.

- ...:

  Not used.

- tol:

  (double) Tolerance parameter used for the checks.

- check_zeroes:

  If `TRUE` then the procedure checks that the given zeroes are zeroes
  of the rational matrix \\x(z)\\.

## Value

Rational matrix object which represents the rational matrix \\x(z)U(z)\\
with the reflected zeroes. The object is of the same class as the input
object `x`.

## Note

The procedures are only implemented for rational matrices with real
coefficients. Therefore complex zeroes occur in complex conjugated pairs
and such pairs are \*\*jointly\*\* reflected to ensure that the result
again has real coefficients. Note however, that the argument `zeroes`
must only contain **one** element of the pairs to be reflected.

In some degnerate cases the procedure(s) may run into numerical
problems: For polynomial matrices and matrices in `lmfd` form, the
procedure assumes that the value of the matrix evaluated at \\z=0\\ is
non singular. For the matrices in statespace form, zeroes on (or close
to) the unit circle are not allowed. In addition multiple zeroes may
lead to unreliable results.

## See also

The procedure uses the helper functions
[`blaschke`](https://bfunovits.github.io/rationalmatrices/reference/blaschke.md)
and
[`make_allpass`](https://bfunovits.github.io/rationalmatrices/reference/make_allpass.md).
For reflecting poles, see
[`reflect_poles`](https://bfunovits.github.io/rationalmatrices/reference/reflect_poles.md).

## Examples

``` r
# ###################################################################
# polynomial matrix 

# construct a (2 x 2) polynomial matrix with degree 2
x = polm(array(c(1,0,
                 0,1,
                 -1,0,
                 0,0,
                 1,0.25,
                 0.75,-1), dim = c(2,2,3)))
x
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2]
#> [1,]        1     0       -1     0     1.00  0.75
#> [2,]        0     1        0     0     0.25 -1.00

# determine zeroes of x(z)
( x_zeroes = zeroes(x) )
#> [1]  0.9236707+0.0000000i  0.4447071-0.8609171i  0.4447071+0.8609171i
#> [4] -0.9709796+0.0000000i

# evaluate the polynomial at the (real) zeroe x_zeroes[1] 
( x_value = zvalue(x, x_zeroes[1]) )
#>              [,1]         [,2]
#> [1,] 0.9294969+0i 0.6398757+0i
#> [2,] 0.2132919+0i 0.1468324+0i

# Verify that the matrix evaluated at z = x_zeroes[1] is rank deficient
d = svd(x_value)$d
(min(d) < (max(d) * sqrt(.Machine$double.eps)))
#> [1] TRUE

# Reflect this zeroe at the unit circle and check the zeroes of the result
x1 = reflect_zeroes(x, x_zeroes[1])

r_zeroes = x_zeroes
r_zeroes[1] = 1 / r_zeroes[1]
x1_zeroes = zeroes(x1)
j = match_vectors(r_zeroes, x1_zeroes)
all.equal(r_zeroes, x1_zeroes[j]) 
#> [1] TRUE

# x_zeroes[2] and x_zeroes[3] form a pair of complex conjugated zeroes
( x_value = zvalue(x, x_zeroes[2]) )
#>                        [,1]                  [,2]
#> [1,]  0.01187902+0.0952052i -0.4075604-0.5742839i
#> [2,] -0.13585347-0.1914280i  1.5434139+0.7657119i

# Verify that the matrix evaluated at x_zeroes[2] is rank deficient
d = svd(x_value)$d
(min(d) < (max(d) * sqrt(.Machine$double.eps)))
#> [1] TRUE

# reflect the real zeroe x_zeroes[1] and this pair of 
# complex conjugated zeroes at the unit circle and verify the result
x1 = reflect_zeroes(x, x_zeroes[c(1,2)])

r_zeroes = x_zeroes
r_zeroes[1:3] = 1 / r_zeroes[1:3]
x1_zeroes = zeroes(x1)
j = match_vectors(r_zeroes, x1_zeroes)
all.equal(r_zeroes, x1_zeroes[j]) 
#> [1] TRUE

# Check that the transformation matrix U (x1 = x %r% U) is all-pass
all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#> [1] TRUE


# ###################################################################
# rational matrix in LMFD form

set.seed(12345)
(x = test_lmfd(dim = c(3,3), degrees = c(2,2), digits = 2))
#> ( 3 x 3 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 2)
#> left factor a(z):
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]        1     0     0     0.59 -0.45  0.63    -0.92  0.37  0.82
#> [2,]        0     1     0     0.71  0.61 -0.28    -0.12  0.52 -0.89
#> [3,]        0     0     1    -0.11 -1.82 -0.28     1.82 -0.75 -0.33
#> right factor b(z):
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]     1.12  1.46 -1.60     0.62  0.81  1.63    -0.32  0.03 -1.06
#> [2,]     0.30 -0.64  1.81     0.61  2.20  0.25    -1.66  1.13  0.94
#> [3,]     0.78 -1.55 -0.48    -0.16  2.05  0.49     1.77 -2.38  0.85
(x_zeroes = zeroes(x))
#> [1] -0.4672627+0.0000000i  0.3378125-0.5740499i  0.3378125+0.5740499i
#> [4]  0.1408149-1.3499406i  0.1408149+1.3499406i  4.9799815+0.0000000i

# reflect all zeroes inside the unit circle
# note: for complex zeroes, select only one of the complex conjugated pair!
x1 = reflect_zeroes(x, x_zeroes[(abs(x_zeroes) <1) & (Im(x_zeroes) >=0 )])

r_zeroes = x_zeroes
r_zeroes[abs(r_zeroes) < 1] = 1 / r_zeroes[abs(r_zeroes) < 1]
(x1_zeroes = zeroes(x1))
#> [1]  0.1408149-1.349941i  0.1408149+1.349941i  0.7614382-1.293923i
#> [4]  0.7614382+1.293923i -2.1401238+0.000000i  4.9799815+0.000000i
j = match_vectors(r_zeroes, x1_zeroes)
all.equal(r_zeroes, x1_zeroes[j]) 
#> [1] TRUE

# Check that the transformation matrix U (x1 = x %r% U) is all-pass
all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#> [1] TRUE

set.seed(NULL)


# ###################################################################
# rational matrix in statespace form (stsp object)

# create a random (2,2) rational matrix in state space form with
# state dimension s=5
set.seed(12345)
(x = test_stsp(dim = c(2,2), s = 5))
#> statespace realization [2,2] with s = 5 states
#>            s[1]       s[2]       s[3]       s[4]       s[5]        u[1]
#> s[1]  0.5855288 -0.2761841 -0.7505320  1.4557851  0.6121235  0.49118828
#> s[2]  0.7094660 -0.2841597  0.8168998 -0.6443284 -0.1623110 -0.32408658
#> s[3] -0.1093033 -0.9193220 -0.8863575 -1.5531374  0.8118732 -1.66205024
#> s[4] -0.4534972 -0.1162478 -0.3315776 -1.5977095  2.1968335  1.76773385
#> s[5]  0.6058875  1.8173120  1.1207127  1.8050975  2.0491903  0.02580105
#> x[1] -1.8179560  0.3706279  0.2987237 -0.4816474  1.6324456  1.00000000
#> x[2]  0.6300986  0.5202165  0.7796219  0.6203798  0.2542712  0.00000000
#>            u[2]
#> s[1]  1.1285108
#> s[2] -2.3803581
#> s[3] -1.0602656
#> s[4]  0.9371405
#> s[5]  0.8544517
#> x[1]  0.0000000
#> x[2]  1.0000000
# zeroes of x(z)
(x_zeroes = zeroes(x))
#> [1]  0.2849467-0.03028325i  0.2849467+0.03028325i -0.3112259-0.40837152i
#> [4] -0.3112259+0.40837152i -0.5438276+0.00000000i

# reflect all unstable zeroes (inside the unit circle)
# note: for complex zeroes, select only one of the complex conjugated pair!
x1 = reflect_zeroes(x, x_zeroes[(abs(x_zeroes) <1) & (Im(x_zeroes) >=0 )])

r_zeroes = x_zeroes
r_zeroes[abs(r_zeroes) < 1] = 1 / r_zeroes[abs(r_zeroes) < 1]
(x1_zeroes = zeroes(x1))
#> [1] -1.838818+0.0000000i -1.180546-1.5490396i -1.180546+1.5490396i
#> [4]  3.470232-0.3688055i  3.470232+0.3688055i
j = match_vectors(r_zeroes, x1_zeroes)
all.equal(r_zeroes, x1_zeroes[j])
#> [1] TRUE

# Check that the transformation matrix U (x1 = x %r% U) is all-pass
all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#> [1] TRUE

set.seed(NULL)
```
