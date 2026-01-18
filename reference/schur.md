# Schur Decomposition

The Schur decomposition of a (real, square) matrix \\A\\ is \$\$A = U S
U'\$\$ where \\U\\ is an orthogonal matrix and `S` is upper
quasi-triangular with 1-by-1 and 2-by-2 blocks on the diagonal. The
2-by-2 blocks correspond to the complex eigenvalues of \\A\\.

## Usage

``` r
schur(A, select = NULL, tol = sqrt(.Machine$double.eps))
```

## Arguments

- A:

  square, non-empty, real matrix.

- select:

  If non NULL, then the Schur decomposition is reordered such that the
  selected cluster of eigenvalues appear in the leading diagonal blocks
  of the quasi-triangular matrix \\S\\. See the details below.

- tol:

  tolerance used to decide whether the target eigenvalues match the
  eigenvalues of \\A\\. See the details below.

## Value

List with components `U`, `S`, `lambda` and `k`, where `lambda` contains
the (computed) eigenvalues of \\A\\ and `k` indicates how many
eigenvalues have been selected (put to the top).

## Details

The optional parameter `select` determines which eigenvalues,
respectively 1-by-1 and 2-by-2 blocks, should be put to the top of the
matrix \\S\\.

- `select='iuc'`: Select the eigenvalues with moduli less than one
  (inside the unit circle).

- `select='ouc'`: Select the eigenvalues with moduli greater than one
  (outside the unit circle).

- `select='lhf'`: Select the eigenvalues with a negative imaginary part
  (left half plane).

- `select='rhf'`: Select the eigenvalues with a positive imaginary part
  (right half plane).

- `select='real'`: Select the real eigenvalues.

- `select='cplx'`: Select the complex eigenvalues.

- `select` is a (complex) vector of eigenvalues. The function checks
  whether the "target eigenvalues", i.e. the entries of `select`, match
  the eigenvalues of the matrix \\A\\. In addition the procedure also
  makes sure that complex eigenvalues are selected in complex conjugated
  pairs.

The function `schur` is simply a wrapper for
[`qz.dgees`](https://rdrr.io/pkg/QZ/man/dd_qz_dgees.html) and
[`qz.dtrsen`](https://rdrr.io/pkg/QZ/man/dd_qz_dtrsen.html).

## Examples

``` r
# generate a "random" 6-by-6 matrix A
m = 6
set.seed(1532)
A = matrix(stats::rnorm(m*m, sd = 0.5), nrow = m, ncol = m)
set.seed(NULL)

# compute the Schur decomposition of A (and its eigenvalues)
out = schur(A)
lambda = out$lambda # eigenvalues of A

# check A = U S U' 
all.equal(A, out$U %*% out$S %*% t(out$U))
#> [1] TRUE
print(out$S)
#>           [,1]        [,2]        [,3]       [,4]         [,5]       [,6]
#> [1,] -1.258119  0.16047666  0.11317140  0.8606362 -0.003078162  0.6654092
#> [2,]  0.000000 -0.07298067 -1.88015298 -0.6360467  0.070634280 -0.3670152
#> [3,]  0.000000  0.43391617 -0.07298067 -0.2158021  0.295439282 -0.2451203
#> [4,]  0.000000  0.00000000  0.00000000  0.7667058  1.205338328  0.7751998
#> [5,]  0.000000  0.00000000  0.00000000 -0.5427939  0.766705822  0.2745867
#> [6,]  0.000000  0.00000000  0.00000000  0.0000000  0.000000000  0.5968981
print(lambda)
#> [1] -1.25811904+0.0000000i -0.07298067+0.9032324i -0.07298067-0.9032324i
#> [4]  0.76670582+0.8088574i  0.76670582-0.8088574i  0.59689812+0.0000000i

# compute an "ordered" Schur decomposition where the eigenvalues 
# inside the unit circle are put to the top of S:
out = schur(A, 'iuc')
print(out$S)
#>             [,1]        [,2]       [,3]       [,4]        [,5]        [,6]
#> [1,] -0.07298067  0.44125325 -0.3372971 -0.1346907  0.17508505 -0.03645365
#> [2,] -1.84889014 -0.07298067 -0.2813143  0.4372502 -0.08687296  0.50003336
#> [3,]  0.00000000  0.00000000  0.5968981 -0.8305954  0.25520954 -0.64478335
#> [4,]  0.00000000  0.00000000  0.0000000 -1.2581190  0.21109806 -0.56034950
#> [5,]  0.00000000  0.00000000  0.0000000  0.0000000  0.76670582  0.45519124
#> [6,]  0.00000000  0.00000000  0.0000000  0.0000000 -1.43730860  0.76670582
print(out$k) # three eigenvalues are inside the unit circle.
#> [1] 3
print(out$lambda) 
#> [1] -0.07298067+0.9032324i -0.07298067-0.9032324i  0.59689812+0.0000000i
#> [4] -1.25811904+0.0000000i  0.76670582+0.8088574i  0.76670582-0.8088574i

# compute an "ordered" Schur decomposition where the eigenvalues 
# lambda[5] and lambda[6] apear in the top. Note that 
# lambda[5] is complex and hence the procedure also selects 
# the conjugate of lambda[5]:  
out = schur(A, lambda[c(6,5)])
print(out$S)
#>            [,1]      [,2]       [,3]       [,4]        [,5]        [,6]
#> [1,]  0.7667058 0.5359624 -0.2917256 -0.4321118 -0.03373457  0.07288711
#> [2,] -1.2207018 0.7667058  1.0416001 -0.6797241 -0.61359022 -0.43178545
#> [3,]  0.0000000 0.0000000  0.5968981 -0.5021335  0.46948321  0.19968389
#> [4,]  0.0000000 0.0000000  0.0000000 -1.2581190 -0.26080155 -0.25951666
#> [5,]  0.0000000 0.0000000  0.0000000  0.0000000 -0.07298067  0.46503579
#> [6,]  0.0000000 0.0000000  0.0000000  0.0000000 -1.75433548 -0.07298067
print(out$k) # three eigenvalues have been selected
#> [1] 3
print(out$lambda) 
#> [1]  0.76670582+0.8088574i  0.76670582-0.8088574i  0.59689812+0.0000000i
#> [4] -1.25811904+0.0000000i -0.07298067+0.9032324i -0.07298067-0.9032324i

if (FALSE) { # \dontrun{
# If the "target" eigenvalues do not match the eigenvalues of A 
# then "schur" throws an error:
out = schur(A, select = lambda[1:3]+ 1)
} # }
```
