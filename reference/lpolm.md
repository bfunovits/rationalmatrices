# Constructor for Laurent Polynomial Matrices

`lpolm` objects represent Laurent polynomial matrices of the form
\$\$a(z) = a\_{-q} z^{-q} + \cdots + a\_{-1} z^{-1} + a_0 + a_1 z +
\cdots + a_p z^p\$\$ If \\a(z)\\ is an \\(m,n)\\ dimensional Laurent
polynomial matrix (i.e. the coefficients \\a_i\\ are \\(m,n)\\
dimensional real or complex valued matrices), then the `lpolm` object
stores the coefficients in an `(m,n,q+p+1)`-dimensional (numeric or
complex) array \\a\_{-q}, \ldots, a\_{-1}, a\_{0}, a\_{1}, \ldots,
a\_{p}\\, together with a class attribute `c("lpolm","ratm")`.  
The constructor function `lpolm(a, min_deg)` takes an integer `min_deg`
(default is zero) and a (numeric or complex) vector, matrix or
3-dimensional array and returns an `lpolm` object.

## Usage

``` r
lpolm(a, min_deg = 0)
```

## Arguments

- a:

  either a (numeric or complex) vector, matrix or 3-D array. A vector is
  coerced to a scalar (i.e. \\(1,1)\\-dimensional) polynomial and a
  matrix gives a Laurent polynomial matrix whose only coefficient matrix
  is of degree `min_deg`.

- min_deg:

  Integer. Default set to zero. Smallest degree in the Laurent
  polynomial. Negative for Laurent polynomials which cannot be coerced
  to
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects
  ([`as.polm`](https://bfunovits.github.io/rationalmatrices/reference/as.polm.md)
  throws an error).

## Value

An object of class `c("lpolm", "ratm")`. When non-negative, the object
is still of this class. However, it can then be coerced to an object of
class `c("polm", "ratm")` using
[`as.polm`](https://bfunovits.github.io/rationalmatrices/reference/as.polm.md)

## Details

Any of the dimensions of the 3-dimensional array may also be zero. In
particular, if the third dimension is zero, then the `lpolm` object is
interpreted as the zero Laurent polynomial.

This class is special in the sense that it is not possible to *upgrade*
it to an
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
or
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
object.

For important methods and functions for this class have a look at the
"see also" section. Note that some functions are only written for
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
objects:
[`degree`](https://bfunovits.github.io/rationalmatrices/reference/degree.md),
[`col_end_matrix`](https://bfunovits.github.io/rationalmatrices/reference/col_end_matrix.md),
normal forms like
[`snf`](https://bfunovits.github.io/rationalmatrices/reference/snf.md),
[`hnf`](https://bfunovits.github.io/rationalmatrices/reference/hnf.md),
[`whf`](https://bfunovits.github.io/rationalmatrices/reference/whf.md).

## See also

- [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

- [`test_lpolm`](https://bfunovits.github.io/rationalmatrices/reference/test_lpolm.md)
  generates random polynomials.

- checks:
  [`is.lpolm`](https://bfunovits.github.io/rationalmatrices/reference/is.md)

- generic S3 methods: [`dim`](https://rdrr.io/r/base/dim.html),
  [`str`](https://rdrr.io/r/utils/str.html) and
  [`print`](https://rdrr.io/r/base/print.html).

- arithmetics:
  [`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md),
  matrix multiplication
  [`%r%`](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md)

- matrix operations:
  [`t.lpolm`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md),
  [`bind`](https://bfunovits.github.io/rationalmatrices/reference/bind.md),
  [`[.lpolm`](https://bfunovits.github.io/rationalmatrices/reference/extract.md),
  `[<-.lpolm`, ...

- [`prune`](https://bfunovits.github.io/rationalmatrices/reference/prune.md)
  discards leading and trailing zero matrices of Laurent polynomials.

- [`get_bwd`](https://bfunovits.github.io/rationalmatrices/reference/get_fwd.md)
  discards all coefficient matrices pertaining to **negative** powers
  and returns a polm object.

- [`get_fwd`](https://bfunovits.github.io/rationalmatrices/reference/get_fwd.md)
  discards all coefficient matrices pertaining to **positive** powers.
  It returns thus an lpolm object.

- [`polm2fwd`](https://bfunovits.github.io/rationalmatrices/reference/polm2fwd.md)
  performs the transformation \\p(z) \rightarrow p(z^{-1})\\ on a polm
  object and returns an lpolm object.

## Examples

``` r
# (1 x 1) Laurent polynomial matrix a(z) =  3z^{-2} + 2z^{-1} + 1
lpolm(3:1, min_deg = -2)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 0, and minimal degree >= -2
#>      z^-2 [,1] z^-1 [,1] z^0 [,1]
#> [1,]         3         2        1

# Non-negative minimal degrees are allowed too (no implicit coercion to polm object)
lpolm(3:1, min_deg = 2)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 4, and minimal degree >= 2
#>      z^2 [,1] z^3 [,1] z^4 [,1]
#> [1,]        3        2        1

lpolm(matrix(1:4,2,2), min_deg = -2)
#> ( 2 x 2 ) Laurent polynomial matrix with degree <= -2, and minimal degree >= -2
#>      z^-2 [,1]  [,2]
#> [1,]         1     3
#> [2,]         2     4
```
