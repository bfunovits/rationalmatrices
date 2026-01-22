# Match Two Vectors

Given two vectors `x,y` of length \\p \leq q\\ respectively, the routine
`match_vectors` returns an integer vector `j`, with unique elements,
such that `x` matches `y[j]` as best as possible. The procedure uses the
"Munkres" algorithm for solving this assignment problem. The procedure
throws an error if the length of `x` is larger than the length of `y`.

## Usage

``` r
match_vectors(x, y = Conj(x))
```

## Arguments

- x, y:

  two vectors of length \\p \leq q\\ respectively.

## Value

The \\p\\-dimensional integer vector `j` (with unique elements).

## Examples

``` r
# Match the roots of two polynomials a1 and a2
p = 5
a1 = rnorm(p+1)
a2 = a1 + rnorm(p+1)*(1e-6) # a2 is a "noisy" copy of a1
z1 = polyroot(a1)
z2 = polyroot(a2)[order(stats::rnorm(p))] # reshuffle the roots of a2
j = match_vectors(z1, z2)
print(data.frame(z1 = z1, j = j, `z2[j]` = z2[j], d = z1-z2[j]))
#>                         z1 j                    z2.j.
#> 1  0.4551862+8.722798e-18i 2  0.4551868+8.723903e-18i
#> 2  0.1034659+9.224606e-01i 3  0.1034656+9.224607e-01i
#> 3  0.1034659-9.224606e-01i 5  0.1034656-9.224607e-01i
#> 4 -1.4531582-1.134652e+00i 1 -1.4531587-1.134653e+00i
#> 5 -1.4531582+1.134652e+00i 4 -1.4531587+1.134653e+00i
#>                             d
#> 1 -6.333816e-07-1.105113e-21i
#> 2  2.543231e-07-1.843963e-07i
#> 3  2.543231e-07+1.843963e-07i
#> 4  4.839707e-07+4.285411e-07i
#> 5  4.839707e-07-4.285411e-07i

# A polynomial with real coefficients has pairs of complex conjugate roots.
# However, the roots returned by "polyroot" in general do not have this 
# property!
# Match the roots and their complex conjugates
j = match_vectors(z1, Conj(z1))
print(data.frame(z = z1, j = j, `Conj(z[j])` = Conj(z1[j]), 
                 d = z1-Conj(z1[j])))
#>                          z j               Conj.z.j..
#> 1  0.4551862+8.722798e-18i 1  0.4551862-8.722798e-18i
#> 2  0.1034659+9.224606e-01i 3  0.1034659+9.224606e-01i
#> 3  0.1034659-9.224606e-01i 2  0.1034659-9.224606e-01i
#> 4 -1.4531582-1.134652e+00i 5 -1.4531582-1.134652e+00i
#> 5 -1.4531582+1.134652e+00i 4 -1.4531582+1.134652e+00i
#>                            d
#> 1  0.000000e+00+1.74456e-17i
#> 2 -5.551115e-17+0.00000e+00i
#> 3  5.551115e-17+0.00000e+00i
#> 4  0.000000e+00+0.00000e+00i
#> 5  0.000000e+00+0.00000e+00i
```
