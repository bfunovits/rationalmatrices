# Purge Rows or Columns of a Polynomial Matrix

This helper function is the main work horse for computing the Hermite
normal form (see
[`hnf`](https://bfunovits.github.io/rationalmatrices/reference/hnf.md))
and the Smith normal form (see
[`snf`](https://bfunovits.github.io/rationalmatrices/reference/snf.md))
of polynomial matrices. It "purges" all elements below, above, to the
right or to the left of a pivot element by elementary row- or column-
operations. Here "purge" means that the elements are either reduced to
zero or that the degree of the elements is made smaller than the degree
of the pivot element.

## Usage

``` r
purge_rc(
  a,
  pivot = c(1, 1),
  direction = c("down", "up", "left", "right"),
  permute = TRUE,
  tol = sqrt(.Machine$double.eps),
  monic = FALSE,
  debug = FALSE
)
```

## Arguments

- a:

  Polynomial matrix of dimension \\(m,n)\\, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

- pivot:

  Integer vector of length 2. Specifies the position of the "pivot"
  element.

- direction:

  Character string.

  down

  :   (default) "purge" all elements below the pivot element by
      elementary row-operations.

  up

  :   "purge" all elements above the pivot element by elementary
      row-operations

  right

  :   "purge" all elements to the right of the pivot element by
      elementary column-operations

  left

  :   "purge" all elements to the left of the pivot element by
      elementary column-operations

- permute:

  Logical, defaults to TRUE. See the details below.

- tol:

  Tolerance parameter, used for "pruning" the polynomial matrix (after
  each step). See
  [`prune`](https://bfunovits.github.io/rationalmatrices/reference/prune.md).

- monic:

  Logical, defaults to FALSE. If TRUE, the coefficient pertaining to the
  highest degree of the pivot element will be normalized to 1.

- debug:

  Logical, default to FALSE. If TRUE, then some diagnostic messages are
  printed.

## Value

List with three slots

- `h`:

  Polynomial matrix of dimension \\(m,n)\\, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).
  This matrix is the result of the "purging" operation(s).

- `u`:

  Unimodular polynomial matrix, i.e. a class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object. For `direction = 'down'` or `direction = 'up'` the matrix
  \\u(z)\\ is \\(m,m)\\ dimensional and satisfies \\a(z) = u(z) h(z)\\.
  For `direction = 'left'` or `direction = 'right'` the matrix \\u(z)\\
  is \\(n,n)\\ dimensional and satisfies \\a(z) = h(z) u(z)\\.

- `u_inv`:

  Unimodular polynomial matrix, i.e. a class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object. This matrix is the inverse of \\u(z)\\.

## Details

Suppose that the matrix \\a(z)\\ has \\m\\ rows, \\n\\ columns and that
the pivot is at position \\(i,j)\\. Furthermore, let us first consider
the case `direction='down'` and `permute=FALSE`. In this case a suitable
multiple - which is computed by the Euclidean polynomial division
algorithm - of the \\i\\-th row is subtracted from all rows below the
\\i\\-th row such that the respective degree of the elements below the
pivot element have a degree which is smaller than the degree of the
pivot.

If the option `permute=TRUE` then first the rows \\i:m\\ are permuted
such that the \\(i,j)\\-th element has the smallest degree among all
elements in the \\j\\-th column and the rows \\i:m\\. Next a suitable
multiple of the \\i\\-th row is subtracted from all rows below the
\\i\\-th row such that the respective degree of the elements below the
pivot element have a degree which is smaller than the degree of the
pivot. These two steps are repeated until all elements below the pivot
element are zero.

Quite analogously the cases `direction='up'`, `direction='left'` and
`direction='right'` may be discussed. Note however, that for the cases
`direction='left'` and `direction='right'` elementary column-operations
are used.

Finally, for `monic=TRUE` the pivot element is made monic, by
multiplying the respective row (or column) by a suitable scalar.

## Examples

``` r
# Generate matrix polynomial
a = test_polm(dim = c(2,3), degree = 1)
print(a, format = 'c')
#> ( 2 x 3 ) matrix polynomial with degree <= 1 
#>             [,1]        [,2]        [,3]
#> [1,]  110 + 111z  120 + 121z  130 + 131z
#> [2,]  210 + 211z  220 + 221z  230 + 231z 

########################################################
# Purge first column downwards
out = purge_rc(a, pivot = c(1,1), monic = TRUE)

# First col zero except for (1,1) element
print(out$h, digits = 2, format = 'c')
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>       [,1]                    [,2]                    [,3]
#> [1,]     1                -9 - 10z               -19 - 20z
#> [2,]     0  1110 + 2220z + 1110z^2  2220 + 4440z + 2220z^2 

# Check polynomial matrix products
all.equal(a, prune(out$u %r% out$h))
#> [1] TRUE
all.equal(out$h, prune(out$u_inv %r% a))
#> [1] TRUE
all.equal(polm(diag(2)), prune(out$u_inv %r% out$u))
#> [1] TRUE

########################################################
# Purge last column upwards
out = purge_rc(a, pivot = c(2,3), direction = "up", monic = TRUE)

# Last col zero except for (2,3) element
print(out$h, digits = 2, format = 'c')
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>                          [,1]                     [,2]  [,3]
#> [1,]  -4620 - 9240z - 4620z^2  -2310 - 4620z - 2310z^2     0
#> [2,]                 21 + 20z                 11 + 10z     1 

# Check polynomial matrix products
all.equal(a, prune(out$u %r% out$h))
#> [1] TRUE
all.equal(out$h, prune(out$u_inv %r% a))
#> [1] TRUE
all.equal(polm(diag(2)), prune(out$u_inv %r% out$u))
#> [1] TRUE

########################################################
# Purge first row right
out = purge_rc(a, pivot = c(1,1), direction = "right", monic = TRUE)

# first row zero except for (1,1) element
print(out$h, digits = 2, format = 'c')
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>             [,1]                       [,2]  [,3]
#> [1,]           1                          0     0
#> [2,]  -99 - 100z  11100 + 22200z + 11100z^2     0 

# Check polynomial matrix products
all.equal(a, prune(out$h %r% out$u))
#> [1] TRUE
all.equal(out$h, prune(a %r% out$u_inv))
#> [1] TRUE
all.equal(polm(diag(3)), prune(out$u_inv %r% out$u))
#> [1] TRUE

########################################################
# Purge last row left
out = purge_rc(a, pivot = c(2,3), direction = "left", monic = TRUE)

# last row zero except for (2,3) element
print(out$h, digits = 2, format = 'c')
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>       [,1]                        [,2]        [,3]
#> [1,]     0  -23100 - 46200z - 23100z^2  101 + 100z
#> [2,]     0                           0           1 

# Check polynomial matrix products
all.equal(a, prune(out$h %r% out$u))
#> [1] TRUE
all.equal(out$h, prune(a %r% out$u_inv))
#> [1] TRUE
all.equal(polm(diag(3)), prune(out$u_inv %r% out$u))
#> [1] TRUE
```
