# Smith Normal Form

Calculates the *Smith normal form* of an \\(m,n)\\-dimensional
polynomial matrix \\a(z)\\, i.e. a factorization of \\a(z)\\ of the form
\$\$a(z) = u(z) s(z) v(z),\$\$ where \\u(z)\\ and \\v(z)\\ are
unimodular polynomial matrices of dimensions \\(m,m)\\ and \\(n,n)\\
respectively, and \\s(z)\\ is a quasi-diagonal polynomial matrix of
dimension \\(m,n)\\ with diagonal elements \\d_i(z)\\ which satisfy

- \\d\_{ii}\\ is monic for \\i\leq r\\ and \\d\_{ii}\\ is zero for
  \\i\>r\\, and

- \\d\_{ii}\\ divides \\d\_{i+1,i+1}\\ for \\i \< r\\.

Here \\0\leq r \leq \min(m,n)\\ is the rank of \\a(z)\\ when considered
as a rational matrix. See \*Gohberg, Lancaster, Rodman 09 - Matrix
Polynomials\* page 318 or (Hannan and Deistler 2012) page 42, Lemma
2.2.3.

## Usage

``` r
snf(a, tol = sqrt(.Machine$double.eps), debug = FALSE)
```

## Arguments

- a:

  Polynomial matrix, i.e. an object of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md).

- tol:

  Tolerance parameter, used for "pruning" the polynomial matrix (after
  each step). See
  [`prune`](https://bfunovits.github.io/rationalmatrices/reference/prune.md).

- debug:

  Logical, default to FALSE. If TRUE, then some diagnostic messages are
  printed.

## Value

A list with five elements

- `s`:

  Matrix polynomial of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  and dimensions \\(m,n)\\, representing \\s(z)\\ in the Smith-form
  (whose only non-zero elements are on the diagonal).

- `u`:

  Unimodular matrix polynomial of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  and dimensions \\(m,m)\\, Represents \\u(z)\\ in the Smith form \\a(z)
  = u(z) s(z) v(z)\\.

- `u_inv`:

  Unimodular matrix polynomial of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  and dimensions \\(m,m)\\, Represents the inverse of \\u(z)\\.

- `v`:

  Unimodular matrix polynomial of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  and dimensions \\(n,n)\\, Represents \\v(z)\\ in the Smith form \\a(z)
  = u(z) s(z) v(z)\\.

- `v_inv`:

  Unimodular matrix polynomial of class
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  and dimensions \\(n,n)\\, Represents the inverse of \\v(z)\\.

## References

Hannan EJ, Deistler M (2012). *The Statistical Theory of Linear
Systems*, Classics in Applied Mathematics. SIAM, Philadelphia.
Originally published: John Wiley & Sons, New York, 1988.

## Examples

``` r
##############
# Quadratic case

a = test_polm(dim = c(2,2), degree = 1)
out = snf(a)

print(out$s, digits = 2, format = 'c')
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>       [,1]          [,2]
#> [1,]     1             0
#> [2,]     0  1 + 2z + z^2 

all.equal(a, prune(out$u %r% out$s %r% out$v))
#> [1] TRUE

##############
# Tall case

a = test_polm(dim = c(3,2), degree = 1)
out = snf(a)

print(out$s, digits = 2, format = 'c')
#> ( 3 x 2 ) matrix polynomial with degree <= 2 
#>       [,1]          [,2]
#> [1,]     1             0
#> [2,]     0  1 + 2z + z^2
#> [3,]     0             0 

all.equal(a, prune(out$u %r% out$s %r% out$v))
#> [1] TRUE

##############
# Wide case

a = test_polm(dim = c(2,3), degree = 1)
out = snf(a)

print(out$s, digits = 2, format = 'c')
#> ( 2 x 3 ) matrix polynomial with degree <= 2 
#>       [,1]          [,2]  [,3]
#> [1,]     1             0     0
#> [2,]     0  1 + 2z + z^2     0 

all.equal(a, prune(out$u %r% out$s %r% out$v))
#> [1] TRUE

##############
# Diagonal case 
z = polm(c(0,1))
a = polm(diag(3))
a[3,3] = 1+2*z 
a[2,2] = a[3,3] * (1-z)
a[1,1] = a[2,2] * (1+z)
print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>                      [,1]          [,2]    [,3]
#> [1,]  1 + 2z - z^2 - 2z^3             0       0
#> [2,]                    0  1 + z - 2z^2       0
#> [3,]                    0             0  1 + 2z 

out = snf(a)

print(out$s, digits = 2, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 3 
#>          [,1]               [,2]                     [,3]
#> [1,]  0.5 + z                  0                        0
#> [2,]        0  -0.5 - 0.5z + z^2                        0
#> [3,]        0                  0  -0.5 - z + 0.5z^2 + z^3 

all.equal(a, prune(out$u %r% out$s %r% out$v))
#> [1] TRUE

##############
# Common factor(s) 
a = test_polm(dim = c(3,3), degree = 1, random = TRUE, digits = 1)
a = a * z
a[,2] = a[,2] * (1+z)
a[,1] = a[,2] * (1-z)

print(a, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 4 
#>                                   [,1]                     [,2]            [,3]
#> [1,]  -0.9z + 3.1z^2 + 0.9z^3 - 3.1z^4  -0.9z + 2.2z^2 + 3.1z^3  -2.1z + 0.1z^2
#> [2,]  -1.1z - 0.4z^2 + 1.1z^3 + 0.4z^4  -1.1z - 1.5z^2 - 0.4z^3           -0.6z
#> [3,]                    -0.2z + 0.2z^3           -0.2z - 0.2z^2   2.5z - 1.2z^2 

out = snf(a)

print(out$s, digits = 2, format = 'c')
#> ( 3 x 3 ) matrix polynomial with degree <= 2 
#>       [,1]     [,2]  [,3]
#> [1,]     z        0     0
#> [2,]     0  z + z^2     0
#> [3,]     0        0     0 

all.equal(a, prune(out$u %r% out$s %r% out$v))
#> [1] TRUE
```
