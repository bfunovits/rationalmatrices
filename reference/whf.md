# Wiener-Hopf Factorization

A (Right-) Wiener-Hopf factorization (R-WHF) of a (square
\\(m,m)\\-dimensional, non singular) polynomial matrix \\A(z)\\ is a
factorization of the form \$\$A(z) = A_f(z) A_0(z) A_b(z) = A_r(z)
A_b(z),\$\$ where

- \\A_b(z)\\:

  is a polynomial matrix which has only zeros outside the unit circle.

- \\A_r(z)\\:

  is a column reduced polynomial matrix with column degrees \\\kappa_i\\
  and all zeroes inside the unit circle.

- \\A_0(z)\\:

  is a diagonal matrix with diagonal entries \\z^{\kappa_i}\\ where
  \\\kappa_i \geq \kappa\_{i+1}\\

- \\A_f(z)\\:

  is a polynomial in \\z^{-1}\\. Note that \\A_f(z) = Ar(z)
  A_0^{-1}(z)\\.

The factors \\A_f(z), A_0(z), A_b(z)\\ are called *forward*, *null* and
*backward* components of \\A(z)\\ and the integers
\\(\kappa_1,\ldots,\kappa_n)\\ are the *partial indices* of \\A(z)\\.  
Similarly, the Left-WHF is defined as \$\$A(z) = A_b(z) A_0(z) A_g(z) =
A_b(z) A_r(z),\$\$ where \\A_r(z)\\ is now row-reduced.  
Note that zeroes on the unit circle are not allowed. In this case the
procedure `whf()` throws an error.  
The Wiener-Hopf factorization plays an important role for the analysis
of linear, rational expectation models. See e.g. (Al-Sadoon 2017) .

## Usage

``` r
whf(a, right_whf = TRUE, tol = sqrt(.Machine$double.eps), debug = FALSE)
```

## Arguments

- a:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object, which represents the polynomial matrix \\A(z)\\.

- right_whf:

  Boolean. Default set to TRUE. If FALSE, then the left WHF \$\$A(z) =
  A_b(z) A_0(z) A_f(z) = A_b(z) A_r(z),\$\$ where the matrix \\A_r(z)\\
  is row-reduced.

- tol:

  Tolerance parameter, used for "pruning" the polynomial matrix (after
  each step). See
  [`prune`](https://bfunovits.github.io/rationalmatrices/reference/prune.md).

- debug:

  Logical, default set to FALSE. If TRUE, then some diagnostic messages
  are printed.

## Value

List with components

- af: A
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  object. Forward component \\A_f\\, a Laurent polynomial whose
  coefficients pertaining to positive powers are zero.

- ab: A
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object. Backward component \\A_b\\

- a0: A
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object. Diagonal matrix with monomials of degrees equal to the partial
  indices \\A_0\\

- ar: A
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object) the column reduced polynomial \\A_r\\

- idx: A vector of integers. partial indices
  \\(\kappa_1,\ldots,\kappa_n)\\

## Details

The algorithm is based on the Smith normal form (SNF),
[snf](https://bfunovits.github.io/rationalmatrices/reference/snf.md),
and a column reduction step, see
[col_reduce](https://bfunovits.github.io/rationalmatrices/reference/col_reduce.md).
An alternative is described in (Gohberg et al. 2003) pages 7ff.

## References

Al-Sadoon MM (2017). “The Linear Systems Approach to Linear Rational
Expectations Models.” *Econometric Theory*, 1-31.
[doi:10.1017/S0266466617000160](https://doi.org/10.1017/S0266466617000160)
. (Gohberg et al. 2003)

## Examples

``` r
set.seed(1234) 

# create test polynomial
a = test_polm(dim = c(3,3), deg = 2, digits = 2, random = TRUE)

# compute WHF and print the result
out = whf(a)
print(out$af, digits = 2, format = 'c')
#> ( 3 x 3 ) Laurent polynomial matrix with degree <= 0, and minimal degree >= -1
#>                   [,1]             [,2]              [,3]
#> [1,]                 0        -2.44z^-1              3.81
#> [2,]  -1.69z^-1 - 2.51  1.97z^-1 + 1.97     -0.42z^-1 - 1
#> [3,]  -0.59z^-1 - 5.78  1.22z^-1 + 5.56  -0.34z^-1 - 1.98 
print(out$a0, digits = 2, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>       [,1]  [,2]  [,3]
#> [1,]     z     0     0
#> [2,]     0     z     0
#> [3,]     0     0     z 
print(out$ab, digits = 2, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>                [,1]          [,2]          [,3]
#> [1,]    1.3 - 5.01z  0.68 + 1.39z  -0.01 + 3.7z
#> [2,]   0.49 - 5.26z  0.96 + 1.48z  0.23 + 3.88z
#> [3,]  -3.61 - 0.22z  0.75 - 0.13z  2.46 - 0.18z 

# check the result
all.equal(a, prune(out$ar %r% out$ab))           # A = Ar * Ab
#> [1] TRUE

# check A(z) = Ab(z^{-1}) A0(z) Ab(z)
# generate random complex z's
z = complex(real = rnorm(10), imaginary = rnorm(10))
a_z  = zvalues(a, z)         # A(z)
ab_z = zvalues(out$ab, z)    # Ab(z)
a0_z = zvalues(out$a0, z)    # A0(z)
af_z = zvalues(out$af, 1/z)  # Af(z^{-1})  
attr(af_z, 'z') = z           # in order to combine the 'zvalues' objects, 
                              # the attribute 'z' must be identical
all.equal(a_z, af_z %r% a0_z %r% ab_z)
#> [1] "Mean relative Mod difference: 7.613902"

all.equal(out$idx, degree(out$ar, 'columns'))    # idx = column degrees of Ar
#> [1] TRUE
all(svd(col_end_matrix(out$ar))$d > 1e-7)     # Ar is column reduced
#> [1] TRUE
abs(zeroes(out$ar, print_message = FALSE))       # Ar has zeroes inside the unit circle
#> [1] 0.3582517 0.4852823 0.4852823
abs(zeroes(out$ab, print_message = FALSE))       # Ab zeroes outside the unit circle
#> [1] 1.934098 1.934098 9.472821

set.seed(NULL)
```
