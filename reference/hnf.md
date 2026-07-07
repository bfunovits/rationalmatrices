# Hermite Normal Form

Calculate the *column Hermite* (default) or *row Hermite form* of a
polynomial matrix \\a(z)\\, by using either (elementary) row operations
(default) or column operations.

## Usage

``` r
hnf(a, from_left = TRUE, tol = sqrt(.Machine$double.eps), debug = FALSE)
```

## Arguments

- a:

  Matrix polynomial, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

- from_left:

  Logical. Default set to TRUE, in which case unimodular
  row-transformations are used to obtain the column Hermite normal form,
  i.e. \\a(z) = u(z) h(z)\\. If FALSE, unimodular column-transformations
  are used to obtain the row Hermite normal form, i.e. \\a(z) = h(z)
  u(z)\\.

- tol:

  Tolerance parameter. Default set to `sqrt(.Machine$double.eps)`.

- debug:

  Logical. If TRUE, some diagnostic messages are printed.

## Value

A list with the following slots.

- `h`:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object which represents the triangular matrix \\h(z)\\. Depending on
  `from_left` the matrix \\h(z)\\ is either quasi-upper- or
  quasi-lower-triangular.

- `u_inv`:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object, which represents the unimodular matrix which transform
  \\a(z)\\ into the desired normal form, i.e. \\h(z) = u^{-1}(z) a(z)\\
  or \\h(z) = a(z)u^{-1}(z)\\.

- `u`:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object, which represents the unimodular matrix \\u(z)\\ such that
  \\a(z) = u(z) h(z)\\ or \\a(z) = h(z)u(z)\\.

## Details

For any \\(m,n)\\ dimensional polynomial matrix \\a(z)\\ with rank \\r\\
(when considered as rational matrix) there exists a unimodular matrix
\\u(z)\\ and indices \\j(1)\<j(2)\<\cdots\<j(r)\\ such that \\h(z)=
u^{-1}(z) a(z)\\ is "quasi-upper-triangular" in the sense that

- \\h\_{ij(i)}\\ is monic (the coefficient pertaining to the highest
  degree is equal to one),

- the elements above \\h\_{ij(i)}\\ have lower polynomial degree than
  \\h\_{ij(i)}\\ and

- \\h\_{i,j}\\ is zero for \\i \>r\\ or \\j \< j(i)\\.

The matrix \\h(z)\\ is called the row Hermite form of \\a(z)\\. The
matrix \\u^{-1}(z)\\ corresponds to the sequence of elementary row
operations which renders \\a(z)\\ into the desired upper-triangular
form.

Quite analogously one may transform the matrix \\a(z)\\ by elementary
column operations into "quasi-lower-triangular" form \\h(z) =
a(z)u^{-1}(z)\\. The corresponding normal form is called *row Hermite
form*.

For a more detailed description, see e.g., (Kailath 1980) (page 375,
Theorem 375) or the package vignette "Rational Matrices".

## References

Kailath T (1980). *Linear Systems*. Prentice Hall, Englewood Cliffs, New
Jersey.

## Examples

``` r
#####################################################################
# Generate polynomial matrix
square = test_polm(dim = c(2,2), degree = 3)
# wide matrix, where all elements have a common factor (1-z)
wide = test_polm(dim = c(2,3), degree = 2) * polm(c(1,-1))
# tall matrix with a "right factor" ((2 x2) random polynomial matrix)
tall = test_polm(dim = c(3,2), degree = 1) %r% 
          test_polm(dim = c(2,2), degree = 1, random = TRUE, digits = 1)

a = tall  # choose one of the above cases

############################
# column Hermite form
out = hnf(a)
print(out$h, digits = 2, format = 'c')
#> ( 3 x 2 ) matrix polynomial with degree <= 4 
#>       [,1]                                 [,2]
#> [1,]     1     2.69 + 3.66z + 2.28z^2 + 1.15z^3
#> [2,]     0  2 + 4.96z + 4.91z^2 + 2.96z^3 + z^4
#> [3,]     0                                    0 

# check result(s)
all.equal(a, prune(out$u %r% out$h))
#> [1] TRUE
all.equal(polm(diag(dim(a)[1])), prune(out$u_inv %r% out$u))
#> [1] TRUE
if (dim(a)[1] == dim(a)[2]) {
  rbind(sort(zeroes(a)), sort(zeroes(out$h)))
}

############################
# row Hermite form
out = hnf(a, from_left = FALSE)
print(out$h, digits = 2, format = 'c')
#> ( 3 x 2 ) matrix polynomial with degree <= 4 
#>                                           [,1]                                  [,2]
#> [1,]                                         1                                     0
#> [2,]     -97.31 - 145.77z - 96.6z^2 - 49.14z^3   2 + 4.96z + 4.91z^2 + 2.96z^3 + z^4
#> [3,]  -195.62 - 291.53z - 193.19z^2 - 98.29z^3  4 + 9.91z + 9.83z^2 + 5.91z^3 + 2z^4 

# check result(s)
all.equal(a, prune(out$h %r% out$u))
#> [1] TRUE
all.equal(polm(diag(dim(a)[2])), prune(out$u_inv %r% out$u))
#> [1] TRUE
if (dim(a)[1] == dim(a)[2]) {
  rbind(sort(zeroes(a)), sort(zeroes(out$h)))
}
```
