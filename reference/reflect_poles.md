# Reflect Poles of Rational Matrices

Given a square, rational matrix \\x(z)\\ and a set of poles of \\x(z)\\,
the function `reflect_poles` constructs an allpass matrix \\U(z)\\ such
that \\x(z)U(z)\\ has poles which are the reciprocals of the selected
poles and the other poles are not changed. Also the zeroes are not
changed, i.e. \\x\\ and \\xU\\ have the same zeroes.

## Usage

``` r
reflect_poles(x, poles, ...)

# S3 method for class 'stsp'
reflect_poles(x, poles, tol = sqrt(.Machine$double.eps), ...)

# S3 method for class 'rmfd'
reflect_poles(
  x,
  poles,
  tol = sqrt(.Machine$double.eps),
  check_poles = TRUE,
  ...
)
```

## Arguments

- x:

  Square, rational matrix object (with real coefficients).

- poles:

  Complex or real vector of poles to be reflected. It is assumed that
  for complex conjugated pairs of poles, only **one** is contained in
  this vector.

- ...:

  Not used.

- tol:

  (double) Tolerance parameter used for the checks.

- check_poles:

  If `TRUE` then the procedure checks that the given poles are poles of
  the rational matrix x(z).

## Value

Rational matrix object which represents the rational matrix \\x(z)U(z)\\
with the reflected poles. The object is of the same class as the input
object `x`.

## Note

The procedures are only implemented for rational matrices with real
coefficients. Therefore complex poles occur in complex conjugated pairs
and such pairs are **jointly** reflected to ensure that the result again
has real coefficients. Note however, that the argument `poles` must only
contain **one** element of the pairs to be reflected.

In some degnerate cases the procedure(s) may run into numerical
problems: For matrices in `rmfd` form, the procedure assumes that the
matrix has no poles at \\z=0\\. For the matrices in statespace form,
poles on (or close to) the unit circle are not allowed. In addition
multiple poles may lead to unreliable results.

## See also

The procedure uses the helper functions
[`blaschke`](https://bfunovits.github.io/rationalmatrices/reference/blaschke.md)
and
[`make_allpass`](https://bfunovits.github.io/rationalmatrices/reference/make_allpass.md).
For reflecting zeroes, see
[`reflect_zeroes`](https://bfunovits.github.io/rationalmatrices/reference/reflect_zeroes.md).

## Examples

``` r
# ###################################################################
# rational matrix in statespace form ('stsp' object)

set.seed(12345)
# create random (2,2) rational matrix in state space form 
# with state dimension s=5
( x = test_stsp(dim = c(2,2), s = 5) )
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
# poles of x(z)
( x_poles = poles(x) )
#> [1]  0.322871350+0.0000000i -0.332609074+0.0000000i  0.009152877-0.8604103i
#> [4]  0.009152877+0.8604103i -4.017541765+0.0000000i

# reflect all unstable poles (inside the unit circle) ###########
# note: for complex zeroes, select only one of the complex conjugated pair!
x1 = reflect_poles(x, poles = x_poles[(abs(x_poles) < 1) & (Im(x_poles) >= 0)])

r_poles = x_poles
r_poles[abs(r_poles) < 1] = 1 / r_poles[abs(r_poles) < 1]
(x1_poles = poles(x1))
#> [1]  0.01236224-1.162105i  0.01236224+1.162105i -3.00653253+0.000000i
#> [4]  3.09720885+0.000000i -4.01754177+0.000000i
j = match_vectors(r_poles, x1_poles)
all.equal(r_poles, x1_poles[j]) 
#> [1] TRUE

# Check that the transformation matrix U (x1 = x %r% U) is all-pass
all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#> [1] TRUE

set.seed(NULL)


# ###################################################################
# rational matrix in RMFD form

set.seed(12345)
# create random (2 x 2) rational matrix in RMFD form with degrees (2,1)
( x = test_rmfd(dim = c(2,2), degree = c(2,1)) )
#> ( 2 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 2, deg(d(z)) = q = 1
#> left factor d(z):
#>        z^0 [,1]       [,2]  z^1 [,1]       [,2]
#> [1,] -0.2841597 -0.1162478 0.3706279 -0.7505320
#> [2,] -0.9193220  1.8173120 0.5202165  0.8168998
#> right factor c(z):
#>      z^0 [,1]  [,2]  z^1 [,1]       [,2]   z^2 [,1]       [,2]
#> [1,]        1     0 0.5855288 -0.1093033  0.6058875  0.6300986
#> [2,]        0     1 0.7094660 -0.4534972 -1.8179560 -0.2761841
# poles of x(z)
( x_poles = poles(x) )
#> [1] -0.492635-0.6478764i -0.492635+0.6478764i  1.045832-0.6704744i
#> [4]  1.045832+0.6704744i

# reflect all unstable poles (inside the unit circle) ###########
# note: for complex zeroes, select only one of the complex conjugated pair!
x1 = reflect_poles(x, poles = x_poles[(abs(x_poles) < 1) & (Im(x_poles) >= 0)])

r_poles = x_poles
r_poles[abs(r_poles) < 1] = 1 / r_poles[abs(r_poles) < 1]
(x1_poles = poles(x1))
#> [1] -0.7436752-0.9780254i -0.7436752+0.9780254i  1.0458317-0.6704744i
#> [4]  1.0458317+0.6704744i
j = match_vectors(r_poles, x1_poles)
all.equal(r_poles, x1_poles[j]) 
#> [1] TRUE

# Check that the transformation matrix U (x1 = x %r% U) is all-pass
all.equal(zvalues(x) %r% Ht(zvalues(x)), zvalues(x1) %r% Ht(zvalues(x1)))
#> [1] TRUE

set.seed(NULL)
```
