# Constructor for Polynomial Matrices

`polm` objects represent polynomial matrices \$\$a(z) = a_0 + a_1 z +
\cdots + a_p z^p\$\$ If the matrix \\a(z)\\ is an \\(m,n)\\ dimensional
polynomial matrix (i.e. the coefficients \\a_i\\ are \\(m,n)\\
dimensional real or complex valued matrices) then the `polm` object
stores the coefficients in an `(m,n,p+1)` dimensional (numeric or
complex) array together with a class attribute `c("polm","ratm")`.  
The constructor function `polm(a)` takes a (numeric or complex) vector,
matrix or 3-dimensional array and returns a `polm` object.

## Usage

``` r
polm(a)
```

## Arguments

- a:

  either a (numeric or complex) vector, matrix or 3-D array. A vector is
  coerced to a scalar (i.e. \\(1,1)\\-dimensional) polynomial and a
  matrix gives a polynomial matrix of zero degree.

## Value

An object of class `polm`.

## Details

Any of the dimensions of the 3-dimensional array may also be zero. In
particular, if the third dimension is zero, then the `polm` object is
interpreted as the zero polynomial.

For important methods and functions for this class have a look at the
"see also" section.

## See also

- [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  objects allow for coefficient matrices pertaining to negative powers
  of \\z\\.

- [`test_polm`](https://bfunovits.github.io/rationalmatrices/reference/test_polm.md)
  generates random polynomials.

- checks:
  [`is.polm`](https://bfunovits.github.io/rationalmatrices/reference/is.md),
  [`is.miniphase`](https://bfunovits.github.io/rationalmatrices/reference/check.md),
  and
  [`is.coprime`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md).
  As a byproduct `is.coprime(a)` computes the zeroes of a square
  polynomial matrix \\a(z)\\.

- generic S3 methods: [`dim`](https://rdrr.io/r/base/dim.html),
  [`str`](https://rdrr.io/r/utils/str.html) and
  [`print`](https://rdrr.io/r/base/print.html).

- arithmetics:
  [`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md),
  matrix multiplication
  [`%r%`](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md),
  polynomial division `%/%`, polynomial remainder `%%`, ...

- matrix operations:
  [`t.polm`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md),
  [`bind`](https://bfunovits.github.io/rationalmatrices/reference/bind.md),
  [`[.polm`](https://bfunovits.github.io/rationalmatrices/reference/extract.md),
  `[<-.polm`, ...

- [`degree`](https://bfunovits.github.io/rationalmatrices/reference/degree.md)
  returns the degree,
  [`col_end_matrix`](https://bfunovits.github.io/rationalmatrices/reference/col_end_matrix.md)
  computes the *column end matrix* and
  [`prune`](https://bfunovits.github.io/rationalmatrices/reference/prune.md)
  "simplifies" a polynomial. Note that the degree of the zero polynomial
  is implemented as being equal to \\-1\\, see
  [`degree`](https://bfunovits.github.io/rationalmatrices/reference/degree.md)!

- [`reflect_zeroes`](https://bfunovits.github.io/rationalmatrices/reference/reflect_zeroes.md)
  may be used to reflect zeroes of a polynomial matrix by multiplication
  with allpass rational matrices.

- normal forms: Hermite normal form
  [`hnf`](https://bfunovits.github.io/rationalmatrices/reference/hnf.md),
  Smith normal form
  [`snf`](https://bfunovits.github.io/rationalmatrices/reference/snf.md),
  column reduced form
  [`col_reduce`](https://bfunovits.github.io/rationalmatrices/reference/col_reduce.md)
  and Wiener Hopf factorization
  [`whf`](https://bfunovits.github.io/rationalmatrices/reference/whf.md).

- [`companion_matrix`](https://bfunovits.github.io/rationalmatrices/reference/companion_matrix.md),
  [`zeroes`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md),
  ...

## Examples

``` r
# (1 x 1) polynomial matrix a(z) =  0 + 1z + 2z^2
polm(0:2)
#> ( 1 x 1 ) matrix polynomial with degree <= 2 
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]        0        1        2

# (2 x 3) polynomial matrix a(z) = a0 (degree is zero)
polm(diag(1, nrow = 2, ncol = 3))
#> ( 2 x 3 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]  [,2]  [,3]
#> [1,]        1     0     0
#> [2,]        0     1     0

# random (2 x 3) polynomial matrix a(z) = a0 + a1 z + a2 z^2 + a3 z^3 (degree = 3)
polm(array(stats::rnorm(2*3*4), dim = c(2,3,4)))
#> ( 2 x 3 ) matrix polynomial with degree <= 3 
#>        z^0 [,1]       [,2]       [,3]  z^1 [,1]      [,2]      [,3]   z^2 [,1]
#> [1,] -0.6775613  2.5608409  0.1927945  1.008282 -1.734546 0.5753295  1.5399118
#> [2,] -0.3507323 -0.1777376 -0.4289415 -0.872292  1.509788 1.0183546 -0.9538422
#>          [,2]      [,3]  z^3 [,1]      [,2]       [,3]
#> [1,] 1.428601 2.0597542 0.1719504 -1.361190 0.54367045
#> [2,] 1.714943 0.9336601 0.5530064  0.768606 0.08115584

# random (2 x 1) polynomial matrix with complex coefficients (degree = 2)
a = polm(array(complex(real = stats::rnorm(2*1*3), 
               imaginary = stats::rnorm(2*1*3)), dim = c(2,1,3)))
is.polm(a)
#> [1] TRUE
dim(a)
#> m n p 
#> 2 1 2 
str(a)
#> ( 2 x 1 ) matrix polynomial with degree <= 2
print(a, digits = 3)
#> ( 2 x 1 ) matrix polynomial with degree <= 2 
#>           z^0 [,1]      z^1 [,1]      z^2 [,1]
#> [1,]  1.279-1.171i -1.440+0.404i -0.136-1.391i
#> [2,] -0.474+1.604i  0.085+0.577i -1.577-1.699i
```
