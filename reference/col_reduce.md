# Construct a Column Reduced Polynomial Matrix

Let \\a(z)\\ be a square (non singular) polynomial matrix. This helper
function constructs a unimodular transformation matrix \\v(z)\\ such
that \\a(z) v^{-1}(z)\\ is column reduced (i.e. the column end matrix
has full rank). Algorithmic implementation are described e.g. in
(Wolovich 1974) Theorem 2.5.7, page 27, (Krishnarao and Chen 1984) , and
(Geurts and Praagman 1996) who show that the KC implementations fails
when the degree of the unimodular matrix \\v^{-1}(z)\\ exceeds the
degree of \\a(z)\\ (page 4 in GP). While all these implementations use
elementary column operations to obtain zero columns in the
column-end-matrix (in order to reduce the degree of the matrix
polynomial), this implementation uses the
[svd](https://rdrr.io/r/base/svd.html). The examples below are taken
from (Geurts and Praagman 1996) and (Krishnarao and Chen 1984) .

## Usage

``` r
col_reduce(a, tol = sqrt(.Machine$double.eps), debug = FALSE)
```

## Arguments

- a:

  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object, which represents the square, polynomial matrix \\a(z)\\.

- tol:

  Double. Tolerance parameter. Default set to
  sqrt(.Machine\$double.eps).

- debug:

  Logical, default to FALSE. If TRUE, then some diagnostic messages are
  printed.

## Value

List with components

- a:

  ([`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  object) is the transformed (column reduced) matrix.

- v,v_inv:

  ([`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects) are the the unimodular matrices \\v(z)\\ and \\v^{-1}(z)\\
  such that \\a(z) v^{-1}(z)\\ is column reduced

.

- col_degrees:

  vector of column degrees of the transformed matrix. Note that the
  columns are permuted such that the transformed matrix has
  \*non-increasing\* column degrees.

- col_end_matrix:

  the column end matrix of the transformed matrix.

## Possible "Improvements"

It is not clear whether the changes are improvements... First: When a
rank deficiency in the column-end-matrix is detected with the SVD, only
one column is "reduced" to zero (even when the rank deficiency is larger
than one). Fewer SVDs are calculated when both It might be better to
reduce all columns to zero which pertain to (numerically) zero
eigenvalues.  
Second: The pivoting mechanism in the QR decomposition might be useful
to single out the columns which should be set to zero. It might be
preferable to the SVD.

## References

Wolovich WA (1974). *Linear multivariable systems*, Applied mathematical
sciences. Springer-Verlag. ISBN 9783540901013. Krishnarao I, Chen CT
(1984). “Two polynomial matrix operations.” *IEEE Transactions on
Automatic Control*, **29**(4), 346-348. Geurts AJ, Praagman C (1996).
“Column Reduction of Polynomial Matrices; Some Remarks on the Algorithm
of Wolovich.” *European Journal of Control*, **2**(2), 152-157.
[doi:10.1016/S0947-3580(96)70039-0](https://doi.org/10.1016/S0947-3580%2896%2970039-0)
.

## See also

The column end matrix may be computed with
[`col_end_matrix`](https://bfunovits.github.io/rationalmatrices/reference/col_end_matrix.md).
The function `col_reduce` is mainly used to compute the Wiener-Hopf
factorization of a polynomial matrix, see
[`whf`](https://bfunovits.github.io/rationalmatrices/reference/whf.md).

## Examples

``` r
# #############################################################################
# define a simple utility function for the computation of the rank of a matrix 
# compare e.g. Matrix::rankMatrix
rkM = function(x) {
  m = nrow(x)
  n = ncol(x)
  tol = max(m,n) * .Machine$double.eps
  
  if (min(m,n) == 0) return(0L)
  s = svd(x, 0, 0)$d
  return(sum(s >= tol*max(s)))
}
# #############################################################################

z = polm(c(0,1))

# Example 2.5.4 in W. A. Wolovich, Linear Multivariable Systems ###############
a = matrix(c(-3,2,0,1,2,3,0,0,2), nrow = 3, ncol = 3) +
    matrix(c(0,4,0,0,0,1,2,0,-3), nrow = 3, ncol = 3) * z + 
    matrix(c(1,0,-1,0,0,0,0,0,0), nrow = 3, ncol = 3) * z^2 
# Original Matrix:
print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]   [,2]    [,3]
#> [1,]  -3 + z^2      1      2z
#> [2,]    2 + 4z      2       0
#> [3,]      -z^2  3 + z  2 - 3z 
# Its column end matrix (and its rank) and column degrees
col_end_matrix(a)
#>      [,1] [,2] [,3]
#> [1,]    1    0    2
#> [2,]    0    0    0
#> [3,]   -1    1   -3
col_end_matrix(a) %>% svd() %>% .$d
#> [1] 3.9516798 0.6198604 0.0000000
degree(a, "c")
#> [1] 2 1 1

# After column reduction:
out = col_reduce(a)
print(out$a, format = 'c', digits = 2)
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>        [,1]    [,2]       [,3]
#> [1,]      1      2z  -3 - 0.5z
#> [2,]      2       0     2 + 3z
#> [3,]  3 + z  2 - 3z      -2.5z 
print(out$col_degrees)
#> [1] 1 1 1
print(out$col_end_matrix)
#>      [,1] [,2] [,3]
#> [1,]    0    2 -0.5
#> [2,]    0    0  3.0
#> [3,]    1   -3 -2.5
print(out$col_end_matrix %>% svd() %>% .$d)
#> [1] 4.6460685 2.7744157 0.4654726

# Check correctness:
all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#> [1] TRUE
all.equal(prune(a), prune(out$a %r% out$v))
#> [1] TRUE


  
# Random example: col degrees = (0,1,-1): throws an error #####################################
a = test_polm(dim = c(3,3), degree = c(0,1,-1), random = TRUE, digits = 2)
print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>        [,1]           [,2]  [,3]
#> [1,]  -2.97   0.54 - 0.55z     0
#> [2,]  -1.57  -1.05 + 0.18z     0
#> [3,]   1.48  -0.62 + 0.83z     0 

if (FALSE) { # \dontrun{
# this throws an error, since a(z) is singular
out = col_reduce(a)
} # }


# Random example: Generic matrices are row reduced ############################################
a = test_polm(dim = c(3,3), degree = c(2,1,0), random = TRUE, digits = 2)
print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                         [,1]           [,2]   [,3]
#> [1,]  1.13 - 0.04z - 0.86z^2  -0.31 + 0.32z  -0.24
#> [2,]  0.56 + 0.86z - 0.35z^2  -0.84 - 1.21z   0.68
#> [3,]  1.79 - 1.16z - 0.59z^2   -0.4 + 0.16z   0.86 

# Column reduction:
out = col_reduce(a) 
print(out$a, format = 'c', digits = 2)
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                         [,1]           [,2]   [,3]
#> [1,]  1.13 - 0.04z - 0.86z^2  -0.31 + 0.32z  -0.24
#> [2,]  0.56 + 0.86z - 0.35z^2  -0.84 - 1.21z   0.68
#> [3,]  1.79 - 1.16z - 0.59z^2   -0.4 + 0.16z   0.86 
print(out$col_degrees)
#> [1] 2 1 0
print(out$col_end_matrix)
#>       [,1]  [,2]  [,3]
#> [1,] -0.86  0.32 -0.24
#> [2,] -0.35 -1.21  0.68
#> [3,] -0.59  0.16  0.86
print(out$col_end_matrix %>% svd() %>% .$d)
#> [1] 1.5309719 1.1339955 0.6573273

# Check:
all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#> [1] TRUE
all.equal(prune(a), prune(out$a %r% out$v))
#> [1] TRUE


# Random example: Column end matrix has rank 2 ################################
col_end_matrix = matrix(round(rnorm(2*3),1), nrow = 3, ncol = 2) %*% 
                 matrix(round(rnorm(2*3),1), nrow = 2, ncol = 3)
a = test_polm(dim = c(3,3), degree = c(2,1,0), random = TRUE, 
               digits = 2, col_end_matrix = col_end_matrix)
print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                         [,1]           [,2]   [,3]
#> [1,]  0.16 - 0.73z - 2.09z^2   1.22 + 0.68z  -0.13
#> [2,]  0.56 - 0.27z + 1.65z^2   0.58 - 0.62z   0.07
#> [3,]  0.08 + 0.71z + 1.69z^2  -1.86 - 0.41z   0.16 
a %>% degree("c")
#> [1] 2 1 0
print(a %>% col_end_matrix)
#>       [,1]  [,2]  [,3]
#> [1,] -2.09  0.68 -0.13
#> [2,]  1.65 -0.62  0.07
#> [3,]  1.69 -0.41  0.16
print(a %>% col_end_matrix %>% svd() %>% .$d)
#> [1] 3.3138389 0.1657469 0.0000000

# Column reduction:
out = col_reduce(a) 
print(out$a, format = 'c', digits = 2)
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>                [,1]          [,2]   [,3]
#> [1,]   1.22 + 0.68z  0.16 + 1.79z  -0.13
#> [2,]   0.58 - 0.62z  0.56 + 0.93z   0.07
#> [3,]  -1.86 - 0.41z  0.08 - 3.13z   0.16 
print(out$col_degrees)
#> [1] 1 1 0
print(out$col_end_matrix)
#>       [,1]       [,2]  [,3]
#> [1,]  0.68  1.7913333 -0.13
#> [2,] -0.62  0.9286667  0.07
#> [3,] -0.41 -3.1340000  0.16
print(out$col_end_matrix %>% svd() %>% .$d)
#> [1] 3.76950806 0.86361997 0.02285663

# Check:
all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#> [1] TRUE
all.equal(prune(a), prune(out$a %r% out$v))
#> [1] TRUE


# Random example: Column end matrix has rank 1 ################################ 
col_end_matrix = matrix(round(rnorm(3),1), nrow = 3, ncol = 1) %*% 
                 matrix(round(rnorm(3),1), nrow = 1, ncol = 3)
a = test_polm(dim = c(3,3), degree = c(2,1,1), random = TRUE, 
               digits = 2, col_end_matrix = col_end_matrix)
print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                          [,1]           [,2]   [,3]
#> [1,]  -0.15 + 0.29z + 0.15z^2   0.17 - 0.03z  -2.48
#> [2,]     1.5 + 0.28z + 0.3z^2   -1.9 - 0.06z  -0.21
#> [3,]   0.61 - 1.62z - 1.35z^2  -1.75 + 0.27z   0.48 
a %>% degree("c")
#> [1] 2 1 0
a %>% col_end_matrix() %>% svd() %>% .$d
#> [1] 2.585712 1.323403 0.000000

# Column reduction:
out = col_reduce(a, debug = FALSE) 
print(out$a, format = 'c', digits = 2)
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>                [,1]           [,2]   [,3]
#> [1,]   0.17 - 0.03z  -0.15 + 1.14z  -2.48
#> [2,]   -1.9 - 0.06z    1.5 - 9.22z  -0.21
#> [3,]  -1.75 + 0.27z  0.61 - 10.37z   0.48 
print(out$col_degrees)
#> [1] 1 1 0
print(out$col_end_matrix)
#>       [,1]   [,2]  [,3]
#> [1,] -0.03   1.14 -2.48
#> [2,] -0.06  -9.22 -0.21
#> [3,]  0.27 -10.37  0.48
print(out$col_end_matrix %>% svd() %>% .$d)
#> [1] 13.9303781  2.4989440  0.2168977

# Check:
all.equal(polm(diag(3)), prune(out$v %r% out$v_inv))
#> [1] TRUE
all.equal(prune(a), prune(out$a %r% out$v))
#> [1] TRUE

#################################################################################

# PG, Ex 1: Jordan Normal Form type (Unimodular matrix) ##########
# Result: Works fine here (does not work with KC implementation)
m0 = diag(3)
m1 = matrix(c(0,1,0,   
              0,0,1,   
              0,0,0), nrow = 3, byrow = TRUE)

# Polymat:
(a = polm(array(c(m0,m1), dim = c(3,3,2))))
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3]
#> [1,]        1     0     0        0     1     0
#> [2,]        0     1     0        0     0     1
#> [3,]        0     0     1        0     0     0
a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>       [,1]  [,2]  [,3]
#> [1,]     1     z     0
#> [2,]     0     1     z
#> [3,]     0     0     1 
a %>% degree("c")
#> [1] 0 1 1

# Column reduced
a_red = col_reduce(a)
a_red$a
#> ( 3 x 3 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]  [,2]  [,3]
#> [1,]        1     0     0
#> [2,]        0     1     0
#> [3,]        0     0     1
degree(a_red$a, "c")
#> [1] 0 0 0

#################################################################################

# PG, Ex 2: Nothing special ##########
# Result: Same as PG

a = test_polm(dim = c(2,2), degree = -1)
a[1,1] = polm(c(4,12,13,6,1))
a[1,2] = polm(c(-2,-5,-4,-1))
a[2,2] = polm(c(2,1))
a
#> ( 2 x 2 ) matrix polynomial with degree <= 4 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2] z^4 [,1]  [,2]
#> [1,]        4    -2       12    -5       13    -4        6    -1        1     0
#> [2,]        0     2        0     1        0     0        0     0        0     0
a %>% print(format = "c")
#> ( 2 x 2 ) matrix polynomial with degree <= 4 
#>                               [,1]                  [,2]
#> [1,]  4 + 12z + 13z^2 + 6z^3 + z^4  -2 - 5z - 4z^2 - z^3
#> [2,]                             0                 2 + z 
col_end_matrix(a)
#>      [,1] [,2]
#> [1,]    1   -1
#> [2,]    0    0
degree(a, "c")
#> [1] 4 3

# Column reduction:
a_red = col_reduce(a)
a_red$a
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]       -2     0       -5     0       -4     0       -1     0
#> [2,]        2     4        1     4        0     1        0     0
a_red$a %>% print(format = "c")
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>                       [,1]          [,2]
#> [1,]  -2 - 5z - 4z^2 - z^3             0
#> [2,]                 2 + z  4 + 4z + z^2 
degree(a_red$a, "c")
#> [1] 3 2

#################################################################################

# PG, Ex 3: Unimodular matrix ##########
# Result: Works, but different matrices, v_inv (U(s) in PG notation is slightly "smaller" here)

a = test_polm(dim = c(3,3), degree = -1)
a[1,1] = polm(c(0,0,0,0,1))
a[1,2] = polm(c(0,0,1))
a[1,3] = polm(c(1,0,0,0,0,0,1))
a[2,1] = polm(c(0,0,1))
a[2,2] = polm(1)
a[2,3] = polm(c(0,0,0,0,1))
a[3,1] = polm(1)
a[3,3] = polm(1)
a
#> ( 3 x 3 ) matrix polynomial with degree <= 6 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3] z^3 [,1]
#> [1,]        0     0     1        0     0     0        0     1     0        0
#> [2,]        0     1     0        0     0     0        1     0     0        0
#> [3,]        1     0     1        0     0     0        0     0     0        0
#>       [,2]  [,3] z^4 [,1]  [,2]  [,3] z^5 [,1]  [,2]  [,3] z^6 [,1]  [,2]  [,3]
#> [1,]     0     0        1     0     0        0     0     0        0     0     1
#> [2,]     0     0        0     0     1        0     0     0        0     0     0
#> [3,]     0     0        0     0     0        0     0     0        0     0     0
a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 6 
#>       [,1]  [,2]     [,3]
#> [1,]   z^4   z^2  1 + z^6
#> [2,]   z^2     1      z^4
#> [3,]     1     0        1 
col_end_matrix(a)
#>      [,1] [,2] [,3]
#> [1,]    1    1    1
#> [2,]    0    0    0
#> [3,]    0    0    0
degree(a, "c")
#> [1] 4 2 6

# Column reduction:
a_red = col_reduce(a)
a_red$a
#> ( 3 x 3 ) matrix polynomial with degree <= 0 
#>      z^0 [,1]  [,2]  [,3]
#> [1,]        0     1     0
#> [2,]        0     0     1
#> [3,]        1     1     0
a_red$a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 0 
#>       [,1]  [,2]  [,3]
#> [1,]     0     1     0
#> [2,]     0     0     1
#> [3,]     1     1     0 
degree(a_red$a, "c")
#> [1] 0 0 0

# Verify the column-reduced matrix time v(z) is equal to original one:
with(a_red, a %r% v) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 6 
#>       [,1]  [,2]     [,3]
#> [1,]   z^4   z^2  1 + z^6
#> [2,]   z^2     1      z^4
#> [3,]     1     0        1 
a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 6 
#>       [,1]  [,2]     [,3]
#> [1,]   z^4   z^2  1 + z^6
#> [2,]   z^2     1      z^4
#> [3,]     1     0        1 

# Unimodular matrix transforming the column-reduced matrix to original one:
a_red$v %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 6 
#>          [,1]  [,2]     [,3]
#> [1,]  1 - z^4  -z^2     -z^6
#> [2,]      z^4   z^2  1 + z^6
#> [3,]      z^2     1      z^4 

# Verify that original matrix a(z) times v^{-1}(z) is column-reduced: PRUNING NECESSARY!
(a %r% a_red$v_inv) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 8 
#>       [,1]                         [,2]                                                [,3]
#> [1,]     0  1 - 8.88178419700125e-16z^6  -2.22044604925031e-16z^2 + 8.88178419700125e-16z^8
#> [2,]     0     -8.88178419700125e-16z^4                         1 + 8.88178419700125e-16z^6
#> [3,]     1                            1                                                   0 
prune(a %r% a_red$v_inv) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 0 
#>       [,1]  [,2]  [,3]
#> [1,]     0     1     0
#> [2,]     0     0     1
#> [3,]     1     1     0 
a_red$v_inv %>% print(format = "c") 
#> ( 3 x 3 ) matrix polynomial with degree <= 6 
#>       [,1]  [,2]           [,3]
#> [1,]     1     0            z^2
#> [2,]  -z^2  -z^4  1 - z^4 + z^6
#> [3,]     0     1           -z^2 

#################################################################################

# PG, Ex 4: PG's "singular case" ##########
# Result: Different results!
#   If eps = 10^(-9) is chosen, 
#   then this algorithm breaks down because the a[0] is recognized as singular! 
#   If eps = 10^(-4) is chosen (as is done in PG), 
#   we obtain a result which is different from the one in PG


# Original matrix (could be argued to be numerically singular at a[0], depending on tolerance!)
a = test_polm(dim = c(3,3), degree = -1)
eps = 10^(-4)
a[1,1] = polm(c(0,0,1,1))
a[1,2] = polm(c(1,eps))
a[1,3] = polm(1)
a[2,1] = polm(c(0,0,2))
a[2,2] = polm(-1)
a[2,3] = polm(-1)
a[3,1] = polm(c(0,0,3))
a[3,2] = polm(1)
a[3,3] = polm(1)
a
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3] z^3 [,1]
#> [1,]        0     1     1        0 1e-04     0        1     0     0        1
#> [2,]        0    -1    -1        0 0e+00     0        2     0     0        0
#> [3,]        0     1     1        0 0e+00     0        3     0     0        0
#>       [,2]  [,3]
#> [1,]     0     0
#> [2,]     0     0
#> [3,]     0     0
a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>            [,1]        [,2]  [,3]
#> [1,]  z^2 + z^3  1 + 1e-04z     1
#> [2,]       2z^2          -1    -1
#> [3,]       3z^2           1     1 
col_end_matrix(a)
#>      [,1]  [,2] [,3]
#> [1,]    1 1e-04    1
#> [2,]    0 0e+00   -1
#> [3,]    0 0e+00    1
col_end_matrix(a) %>% svd() %>% .$d
#> [1] 1.847759e+00 7.653669e-01 1.355253e-20
degree(a, "c")
#> [1] 3 1 0

# Column reduction: Note that there is a column with "small length". 
#   It depends on the tolerance whether this column is considered to be zero.
#   Also, note the singular values of the column-end-matrix of the reduced polymat!     
a_red = col_reduce(a)
a_red$a
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]        0     1     1        0 1e-04     0    -9999     0     0
#> [2,]        0    -1    -1        0 0e+00     0    10002     0     0
#> [3,]        0     1     1        0 0e+00     0    -9997     0     0
a_red$a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                       [,1]        [,2]  [,3]
#> [1,]  -9999.00000001665z^2  1 + 1e-04z     1
#> [2,]   10002.0000000166z^2          -1    -1
#> [3,]  -9997.00000001665z^2           1     1 
a_red$a %>% print(format = "c", digits = 3)
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]  [,2]  [,3]
#> [1,]  -9999z^2     1     1
#> [2,]  10002z^2    -1    -1
#> [3,]  -9997z^2     1     1 

a_red$a %>% col_end_matrix()
#>       [,1]  [,2] [,3]
#> [1,] -9999 1e-04    1
#> [2,] 10002 0e+00   -1
#> [3,] -9997 0e+00    1
a_red$a %>% col_end_matrix() %>% svd() %>% .$d
#> [1] 1.731935e+04 3.560566e-04 8.108103e-05
a_red$a %>% degree("c")
#> [1] 2 1 0

# Col-reduced matrix time v(z) = original: 
# It works up to a small numerical issue in the (1,1) element 
with(a_red, a %r% v) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>                            [,1]        [,2]  [,3]
#> [1,]  z^2 + 1.00000000000166z^3  1 + 1e-04z     1
#> [2,]                       2z^2          -1    -1
#> [3,]                       3z^2           1     1 
a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>            [,1]        [,2]  [,3]
#> [1,]  z^2 + z^3  1 + 1e-04z     1
#> [2,]       2z^2          -1    -1
#> [3,]       3z^2           1     1 

# Check unimodular matrix taking the col-reduced matrix back to original:
# Small (i.e. unproblematic) numerical issue in the (2,1)-element) which is also reflected above.
a_red$v %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                      [,1]  [,2]  [,3]
#> [1,]                    1     0     0
#> [2,]  10000.0000000166z^2     1     0
#> [3,]                    0     0     1 

# Original times v^{-1}(z) = col-reduced:
# Works fine, up to small (non-problematic) issue in (1,3)-element
a_red$a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                       [,1]        [,2]  [,3]
#> [1,]  -9999.00000001665z^2  1 + 1e-04z     1
#> [2,]   10002.0000000166z^2          -1    -1
#> [3,]  -9997.00000001665z^2           1     1 
(a %r% a_red$v_inv) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>                                                 [,1]        [,2]  [,3]
#> [1,]  -9999.00000001665z^2 - 1.66489044772788e-12z^3  1 + 1e-04z     1
#> [2,]                             10002.0000000166z^2          -1    -1
#> [3,]                            -9997.00000001665z^2           1     1 

# Matrix transforming the original to column-reduced:
# Small (unproblematic) issue in (2,3)-element
a_red$v_inv %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                       [,1]  [,2]  [,3]
#> [1,]                     1     0     0
#> [2,]  -10000.0000000166z^2     1     0
#> [3,]                     0     0     1 

# Algebraic result given in PG: Different from the one obtained here! ####
# The element (3,1) of matrix R(s) in PG seems to be incorrect. 
# Changing this element below results in a column-reduced matrix.
eta = 1/eps

r = test_polm(dim = c(3,3), degree = -1)
r[1,1] = polm(c(0,-3*eta))
r[1,2] = polm(c(1))
r[1,3] = polm(1)
r[2,1] = polm(c(0,3*eta))
r[2,2] = polm(c(-1, eps))
r[2,3] = polm(-1)
r[3,1] = polm(c(1,-3*eta,5))
r[3,2] = polm(c(1,-eps))
r[3,3] = polm(1)
r
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]   [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]        0     1     1   -30000  0e+00     0        0     0     0
#> [2,]        0    -1    -1    30000  1e-04     0        0     0     0
#> [3,]        1     1     1   -30000 -1e-04     0        5     0     0
# PG: Col-reduced
r %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                    [,1]         [,2]  [,3]
#> [1,]            -30000z            1     1
#> [2,]             30000z  -1 + 1e-04z    -1
#> [3,]  1 - 30000z + 5z^2   1 - 1e-04z     1 
# Here: Col-reduced
a_red$a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                       [,1]        [,2]  [,3]
#> [1,]  -9999.00000001665z^2  1 + 1e-04z     1
#> [2,]   10002.0000000166z^2          -1    -1
#> [3,]  -9997.00000001665z^2           1     1 

u = test_polm(dim = c(3,3), degree = -1)
u[1,1] = polm(c(1))
u[1,2] = polm(c(0))
u[1,3] = polm(0)
u[2,1] = polm(c(0,-3*eta,-eta))
u[2,2] = polm(c(1))
u[2,3] = polm(0)
u[3,1] = polm(c(0,0,eta+2))
u[3,2] = polm(c(0,-eps))
u[3,3] = polm(1)
u
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]   [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]        1     0     0        0  0e+00     0        0     0     0
#> [2,]        0     1     0   -30000  0e+00     0   -10000     0     0
#> [3,]        0     0     1        0 -1e-04     0    10002     0     0
u %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                     [,1]     [,2]  [,3]
#> [1,]                   1        0     0
#> [2,]  -30000z - 10000z^2        1     0
#> [3,]            10002z^2  -1e-04z     1 

# Check whether their result makes sense:
# Column-reduced matrix r is different from a %r% u!!! 
# This mistake can be corrected by adjusting element (3,1) of r
a %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>            [,1]        [,2]  [,3]
#> [1,]  z^2 + z^3  1 + 1e-04z     1
#> [2,]       2z^2          -1    -1
#> [3,]       3z^2           1     1 

a %r% u
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]   [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]        0     1     1   -30000  0e+00     0        0     0     0
#> [2,]        0    -1    -1    30000  1e-04     0        0     0     0
#> [3,]        0     1     1   -30000 -1e-04     0        5     0     0
(a %r% u) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                 [,1]         [,2]  [,3]
#> [1,]         -30000z            1     1
#> [2,]          30000z  -1 + 1e-04z    -1
#> [3,]  -30000z + 5z^2   1 - 1e-04z     1 
(a %r% u) %>% col_end_matrix()
#>      [,1]   [,2] [,3]
#> [1,]    0  0e+00    1
#> [2,]    0  1e-04   -1
#> [3,]    5 -1e-04    1
(a %r% u) %>% col_end_matrix() %>% svd() %>% .$d
#> [1] 5.107156e+00 1.384541e+00 7.071068e-05
(a %r% u) %>% degree("c")
#> [1] 2 1 0

r
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]   [,2]  [,3] z^2 [,1]  [,2]  [,3]
#> [1,]        0     1     1   -30000  0e+00     0        0     0     0
#> [2,]        0    -1    -1    30000  1e-04     0        0     0     0
#> [3,]        1     1     1   -30000 -1e-04     0        5     0     0
r %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                    [,1]         [,2]  [,3]
#> [1,]            -30000z            1     1
#> [2,]             30000z  -1 + 1e-04z    -1
#> [3,]  1 - 30000z + 5z^2   1 - 1e-04z     1 

#########################
# Change an element in r
r2 = r
r2[3,1] = polm(c(0,-3*eta,5))
r2 %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                 [,1]         [,2]  [,3]
#> [1,]         -30000z            1     1
#> [2,]          30000z  -1 + 1e-04z    -1
#> [3,]  -30000z + 5z^2   1 - 1e-04z     1 
r2 %>% col_end_matrix()
#>      [,1]   [,2] [,3]
#> [1,]    0  0e+00    1
#> [2,]    0  1e-04   -1
#> [3,]    5 -1e-04    1
r2 %>% col_end_matrix() %>% svd() %>% .$d
#> [1] 5.107156e+00 1.384541e+00 7.071068e-05
r2 %>% col_end_matrix() %>% rkM()
#> [1] 3
r2 %>% degree("c")
#> [1] 2 1 0

(a %r% u) %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>                 [,1]         [,2]  [,3]
#> [1,]         -30000z            1     1
#> [2,]          30000z  -1 + 1e-04z    -1
#> [3,]  -30000z + 5z^2   1 - 1e-04z     1 
 

#################################################################################

# PG, Ex 5: Breaks down in PG, but they give an algebraic solution ##########
# Result: Works here! 
#  Different result as algebraic solution indicated in PG (up to column permutation).
#  Of course, the obtained result is also column-reduced!
#  The unimodular matrix (which column-reduces the original polynomial matrix) is different 

a = test_polm(dim = c(4,4), degree = -1)
eps = 10^(-9)
a[1,1] = polm(c(1,1,1))
a[1,2] = polm(c(0,eps))
a[1,3] = polm(c(0,0,0,1))
a[1,4] = polm(c(1,0,0,1))
a[2,1] = polm(c(0,2))
a[2,2] = polm(0)
a[2,3] = polm(c(1,0,0,2))
a[2,4] = a[2,3]
a[3,1] = polm(c(1,3))
a[3,2] = polm(3)
a[3,3] = polm(c(3,0,3))
a[3,4] = a[3,3]
a[4,1] = polm(c(0,4))
a[4,2] = polm(c(0))
a[4,3] = polm(c(1,0,0,4))
a[4,4] = a[4,3]
a
#> ( 4 x 4 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2]  [,3]  [,4] z^1 [,1]  [,2]  [,3]  [,4] z^2 [,1]  [,2]  [,3]
#> [1,]        1     0     0     1        1 1e-09     0     0        1     0     0
#> [2,]        0     0     1     1        2 0e+00     0     0        0     0     0
#> [3,]        1     3     3     3        3 0e+00     0     0        0     0     3
#> [4,]        0     0     1     1        4 0e+00     0     0        0     0     0
#>       [,4] z^3 [,1]  [,2]  [,3]  [,4]
#> [1,]     0        0     0     1     1
#> [2,]     0        0     0     2     2
#> [3,]     3        0     0     0     0
#> [4,]     0        0     0     4     4

# Original matrix, its column-end-matrix with its singular values, and its column degrees
a %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 3 
#>              [,1]    [,2]      [,3]      [,4]
#> [1,]  1 + z + z^2  1e-09z       z^3   1 + z^3
#> [2,]           2z       0  1 + 2z^3  1 + 2z^3
#> [3,]       1 + 3z       3  3 + 3z^2  3 + 3z^2
#> [4,]           4z       0  1 + 4z^3  1 + 4z^3 
col_end_matrix(a)
#>      [,1]  [,2] [,3] [,4]
#> [1,]    1 1e-09    1    1
#> [2,]    0 0e+00    2    2
#> [3,]    0 0e+00    0    0
#> [4,]    0 0e+00    4    4
col_end_matrix(a) %>% svd() %>% .$d
#> [1] 6.484499e+00 9.753345e-01 5.849050e-25 0.000000e+00
col_end_matrix(a) %>% rkM()
#> [1] 2
degree(a, "c")
#> [1] 2 1 3 3

# Column-reduction:
a_red = col_reduce(a)
a_red$a
#> ( 4 x 4 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2]  [,3]  [,4] z^1 [,1]  [,2]  [,3]  [,4]
#> [1,]        1     0     1     0        1     0     0     0
#> [2,]        0     0     0     1        2     0     0     0
#> [3,]        1     3     0     3        3     0     0     0
#> [4,]        0     0     0     1        4     0     0     0
a_red$a %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 1 
#>         [,1]  [,2]  [,3]  [,4]
#> [1,]   1 + z     0     1     0
#> [2,]      2z     0     0     1
#> [3,]  1 + 3z     3     0     3
#> [4,]      4z     0     0     1 
col_end_matrix(a_red$a)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    1    0
#> [2,]    2    0    0    1
#> [3,]    3    3    0    3
#> [4,]    4    0    0    1
col_end_matrix(a_red$a) %>% svd() %>% .$d
#> [1] 6.5378222 2.6726142 0.9977857 0.3441475
degree(a_red$a, "c")
#> [1] 1 0 0 0

# Check whether col-reduced matrix times unimodular v(z) = original:
# Works with small (non-problematic) numerical mistakes.
# Brutal pruning "solves" it
# Similar for the unimodular matrix v(z) itself
with(a_red, a %r% v) %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 3 
#>              [,1]  [,2]                                                   [,3]                                                       [,4]
#> [1,]  1 + z + z^2     0  1.73580252387356e-16z - 4.84643525376756e-17z^2 + z^3  1 + 1.73580252387356e-16z - 4.84643525376756e-17z^2 + z^3
#> [2,]           2z     0                     1 + 3.47160504774711e-16z^2 + 2z^3                         1 + 3.47160504774711e-16z^2 + 2z^3
#> [3,]       1 + 3z     3                       3 + 1.73580252387356e-16z + 3z^2                           3 + 1.73580252387356e-16z + 3z^2
#> [4,]           4z     0                     1 + 6.94321009549423e-16z^2 + 4z^3                         1 + 6.94321009549423e-16z^2 + 4z^3 
with(a_red, a %r% v) %>% prune(brutal = TRUE) %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 3 
#>              [,1]  [,2]      [,3]      [,4]
#> [1,]  1 + z + z^2     0       z^3   1 + z^3
#> [2,]           2z     0  1 + 2z^3  1 + 2z^3
#> [3,]       1 + 3z     3  3 + 3z^2  3 + 3z^2
#> [4,]           4z     0  1 + 4z^3  1 + 4z^3 
a_red$v %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 3 
#>       [,1]  [,2]                            [,3]                               [,4]
#> [1,]     1     0     1.73580252387356e-16z + z^2        1.73580252387356e-16z + z^2
#> [2,]     0     1      0.666666666666666z^2 - z^3         0.666666666666666z^2 - z^3
#> [3,]   z^2     0  -z^2 + 1.92296268638356e-16z^3  1 - z^2 + 1.92296268638356e-16z^3
#> [4,]     0     0                               1                                  1 
a_red$v %>% prune(brutal = TRUE) %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 3 
#>       [,1]  [,2]                        [,3]                        [,4]
#> [1,]     1     0                         z^2                         z^2
#> [2,]     0     1  0.666666666666666z^2 - z^3  0.666666666666666z^2 - z^3
#> [3,]   z^2     0                        -z^2                     1 - z^2
#> [4,]     0     0                           1                           1 

# Check whether original times v^{-1}(z) = col-reduced:
# Same as above: non-problematic numerical mistakes, "solved" by brutally pruning
# Same for unimodular v^{-1}(z)
(a %r% a_red$v_inv) %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 4 
#>                                                             [,1]    [,2]  [,3]                                                                                                  [,4]
#> [1,]                             1 + z - 2.22044604925031e-16z^2  1e-09z     1  -1.73580252387356e-16z + 4.84643525376756e-17z^2 - 6.66667165916124e-10z^3 + 1.00000030478498e-09z^4
#> [2,]                                                          2z       0     0                                                 1 - 3.47160504774711e-16z^2 - 4.44089209850063e-16z^3
#> [3,]  1 + 3z + 2.22044604925031e-16z^2 + 2.22044604925031e-16z^4       3     0                                                   3 - 1.73580252387356e-16z + 8.88178419700125e-16z^2
#> [4,]                                                          4z       0     0                                                 1 - 6.94321009549423e-16z^2 - 8.88178419700125e-16z^3 
(a %r% a_red$v_inv) %>% prune(brutal = TRUE) %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 1 
#>         [,1]  [,2]  [,3]  [,4]
#> [1,]   1 + z     0     1     0
#> [2,]      2z     0     0     1
#> [3,]  1 + 3z     3     0     3
#> [4,]      4z     0     0     1 
a_red$a %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 1 
#>         [,1]  [,2]  [,3]  [,4]
#> [1,]   1 + z     0     1     0
#> [2,]      2z     0     0     1
#> [3,]  1 + 3z     3     0     3
#> [4,]      4z     0     0     1 
a_red$v_inv %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 4 
#>       [,1]  [,2]  [,3]                                     [,4]
#> [1,]     1     0     0             -1.73580252387356e-16z - z^2
#> [2,]     0     1     0              -0.666666666666666z^2 + z^3
#> [3,]   z^2     0    -1  1 - z^2 + 1.87160162510007e-17z^3 - z^4
#> [4,]  -z^2     0     1      z^2 - 1.87160162510007e-17z^3 + z^4 
a_red$v_inv %>% prune(brutal = TRUE) %>% print(format = "c")
#> ( 4 x 4 ) matrix polynomial with degree <= 4 
#>       [,1]  [,2]  [,3]                         [,4]
#> [1,]     1     0     0                         -z^2
#> [2,]     0     1     0  -0.666666666666666z^2 + z^3
#> [3,]   z^2     0    -1                1 - z^2 - z^4
#> [4,]  -z^2     0     1                    z^2 + z^4 

# Krishnarao and Chen example #####################

(a = polm(array(c(4,0,-2,2,  
                  12, 0, -5, 1,    
                  13,0,1,0,   
                  2,0,0,0), dim = c(2,2,4))))
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]        4    -2       12    -5       13     1        2     0
#> [2,]        0     2        0     1        0     0        0     0
out = col_reduce(a)
out$a
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2]
#> [1,]       -2     4       -5    16        1    23
#> [2,]        2     0        1    -4        0    -2
a %r% out$v_inv
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2]
#> [1,]       -2     4       -5    16        1    23
#> [2,]        2     0        1    -4        0    -2
with(out, a %r% v)
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2] z^3 [,1]  [,2]
#> [1,]        4    -2       12    -5       13     1        2     0
#> [2,]        0     2        0     1        0     0        0     0

# Majid: Ex 3 ################

(m0 = matrix(c(-3,1,0, 2,2,0,   0,3,2), nrow = 3, byrow = TRUE))
#>      [,1] [,2] [,3]
#> [1,]   -3    1    0
#> [2,]    2    2    0
#> [3,]    0    3    2
(m1 = matrix(c(0,0,2,   4,0,0,   0,1,-3), nrow = 3, byrow = TRUE))
#>      [,1] [,2] [,3]
#> [1,]    0    0    2
#> [2,]    4    0    0
#> [3,]    0    1   -3
(m2 = matrix(c(1,0,0,   0,0,0,   -1,0,0), nrow = 3, byrow = TRUE))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    0    0
#> [3,]   -1    0    0
m = polm(array(c(m0,m1,m2), dim = c(3,3,3)))
m %>% print(format = "c")
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>           [,1]   [,2]    [,3]
#> [1,]  -3 + z^2      1      2z
#> [2,]    2 + 4z      2       0
#> [3,]      -z^2  3 + z  2 - 3z 
m %>% col_end_matrix()
#>      [,1] [,2] [,3]
#> [1,]    1    0    2
#> [2,]    0    0    0
#> [3,]   -1    1   -3
m %>% col_end_matrix() %>% svd() %>% .$d
#> [1] 3.9516798 0.6198604 0.0000000
m %>% degree("c")
#> [1] 2 1 1

(out = col_reduce(m))
#> $a
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3]
#> [1,]        1     0    -3        0     2  -0.5
#> [2,]        2     0     2        0     0   3.0
#> [3,]        3     2     0        1    -3  -2.5
#> 
#> $v
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3]
#> [1,]        0     1     0      0.5     0     0
#> [2,]        0     0     1      0.5     0     0
#> [3,]        1     0     0      0.0     0     0
#> 
#> $v_inv
#> ( 3 x 3 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2]  [,3] z^1 [,1]  [,2]  [,3]
#> [1,]        0     0     1        0     0   0.0
#> [2,]        1     0     0        0     0  -0.5
#> [3,]        0     1     0        0     0  -0.5
#> 
#> $col_degrees
#> [1] 1 1 1
#> 
#> $col_end_matrix
#>      [,1] [,2] [,3]
#> [1,]    0    2 -0.5
#> [2,]    0    0  3.0
#> [3,]    1   -3 -2.5
#> 

 # Majid: Ex 4 ################3
(m0 = matrix(c(1,0,   0,1), nrow = 2, byrow = TRUE))
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
(m1 = matrix(c(0,0,   2,0), nrow = 2, byrow = TRUE))
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    2    0
(m2 = matrix(c(1,1,   0,0), nrow = 2, byrow = TRUE))
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    0    0
(m3 = matrix(c(2,0,   0,0), nrow = 2, byrow = TRUE))
#>      [,1] [,2]
#> [1,]    2    0
#> [2,]    0    0
m = polm(array(c(m0,m1,m2,m3), dim = c(2,2,4)))
m %>% print(format = "c")
#> ( 2 x 2 ) matrix polynomial with degree <= 3 
#>                 [,1]  [,2]
#> [1,]  1 + z^2 + 2z^3   z^2
#> [2,]              2z     1 
m %>% col_end_matrix()
#>      [,1] [,2]
#> [1,]    2    1
#> [2,]    0    0
m %>% col_end_matrix() %>% svd() %>% .$d
#> [1] 2.236068 0.000000
m %>% degree("c")
#> [1] 3 2

(out = col_reduce(m))
#> $a
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2] z^2 [,1]  [,2]
#> [1,]        0     1        0     0        1     0
#> [2,]        1    -1        0     0        0     0
#> 
#> $v
#> ( 2 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        1     1        2     0
#> [2,]        1     0        0     0
#> 
#> $v_inv
#> ( 2 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]        0     1        0     0
#> [2,]        1    -1        0    -2
#> 
#> $col_degrees
#> [1] 2 0
#> 
#> $col_end_matrix
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    0   -1
#> 
```
