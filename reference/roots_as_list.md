# Polynomial Roots as List

This helper function coerces a vector of (complex or real) roots into a
list. For a complex conjugate pair of roots, only the one with a
positive imaginary part is retained.

## Usage

``` r
roots_as_list(roots, tol = sqrt(.Machine$double.eps))
```

## Arguments

- roots:

  Vector of cplx or doubles. Obtained e.g. from a call to
  [`zeroes`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md).

- tol:

  Double. Tolerance parameter used to decide whether roots and conjugate
  roots "match".

## Value

List of roots. For complex conjugated pairs of roots, only the ones with
positive part are returned. The procedure orders the roots according to
their imaginary part, and thus real roots come first.

## Details

The routine assumes that complex roots appear, up to some small
numerical errors, in complex conjugate pairs. Therefore the procedure
tries to match the roots with their complex conjugates. If this is not
possible, then an error is thrown.

Roots, which are classified as real, are replaced by their real part.
Roots, which are classified as complex, are replaced by the mean of the
root and the best matching conjugate root.

## See also

The matching of conjugate roots is done by the internal helper function
[`match_vectors`](https://bfunovits.github.io/rationalmatrices/reference/match_vectors.md).

## Examples

``` r
set.seed(12345)
p = 5
a = rnorm(p+1)   # coefficients of a random polynomial a(z) of degree p = 5
z = polyroot(a)  # compute the roots of a(z)
z
#> [1]  0.2353419+9.290992e-01i -0.5313699+3.107316e-01i -0.5313699-3.107316e-01i
#> [4]  0.9253355-8.751787e-16i  0.2353419-9.290992e-01i

# try to match roots and conjugate roots
j = match_vectors(z, Conj(z))
# z is approximately equal to Conj( z[j] )
print(data.frame(z = z, j = j, `Conj(z[j])` = Conj(z[j]), d = z-Conj(z[j])))
#>                          z j               Conj.z.j..
#> 1  0.2353419+9.290992e-01i 5  0.2353419+9.290992e-01i
#> 2 -0.5313699+3.107316e-01i 3 -0.5313699+3.107316e-01i
#> 3 -0.5313699-3.107316e-01i 2 -0.5313699-3.107316e-01i
#> 4  0.9253355-8.751787e-16i 4  0.9253355+8.751787e-16i
#> 5  0.2353419-9.290992e-01i 1  0.2353419-9.290992e-01i
#>                             d
#> 1  1.221245e-15-9.992007e-16i
#> 2  3.108624e-15+1.887379e-15i
#> 3 -3.108624e-15+1.887379e-15i
#> 4  0.000000e+00-1.750357e-15i
#> 5 -1.221245e-15-9.992007e-16i
# z[1] and z[5] are complex conjugates (up to numerical errors)
# z[2] and z[3] are complex conjugates (up to numerical errors)
# z[4] is real (up to numerical errors)

# coerce the vector "z" to a list
(z_list = roots_as_list(z))
#> [[1]]
#> [1] 0.9253355
#> 
#> [[2]]
#> [1] -0.5313699+0.3107316i
#> 
#> [[3]]
#> [1] 0.2353419+0.9290992i
#> 
# the first slot contains the real root (and thus is of class "numeric")
# z_list[[1]] = Re(z[4])
# the slots 2 and 3 contain the complex roots and thus are of class "complex"
# z_list[[2]] = (z[2] + Conj(z[3])/2
# z_list[[3]] = (z[1] + Conj(z[5])/2

# The routine zeroes() uses the function eigen() 
# (to compute the eigenvalues of the companion matrix) 
# and thus returns exact conjugate pairs:
(z = zeroes(polm(a)))
#> [1] -0.5313699-0.3107316i -0.5313699+0.3107316i  0.9253355+0.0000000i
#> [4]  0.2353419-0.9290992i  0.2353419+0.9290992i

# match roots and conjugate roots
j = match_vectors(z, Conj(z))
# z is equal to Conj(z[j])
print(j)
#> [1] 2 1 3 5 4
all.equal(z, Conj(z[j]))
#> [1] TRUE

# coerce the vector "z" to a list
(z_list = roots_as_list(z))
#> [[1]]
#> [1] 0.9253355
#> 
#> [[2]]
#> [1] -0.5313699+0.3107316i
#> 
#> [[3]]
#> [1] 0.2353419+0.9290992i
#> 

set.seed(NULL)
```
