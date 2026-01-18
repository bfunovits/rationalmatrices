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
#> 1  0.1081845-3.805031e-23i 1  0.1081838-3.805031e-23i
#> 2  0.6907488+7.141664e-01i 3  0.6907488+7.141650e-01i
#> 3  0.6907488-7.141664e-01i 2  0.6907488-7.141650e-01i
#> 4 -1.9146323+7.502615e-15i 4 -1.9146367+7.530065e-15i
#> 5  7.3319114+1.630502e-15i 5  7.3321159+1.735672e-15i
#>                             d
#> 1  7.160299e-07+0.000000e+00i
#> 2  1.598535e-08+1.405276e-06i
#> 3  1.598535e-08-1.405276e-06i
#> 4  4.372101e-06-2.745073e-17i
#> 5 -2.044562e-04-1.051698e-16i

# A polynomial with real coefficients has pairs of complex conjugate roots.
# However, the roots returned by "polyroot" in general do not have this 
# property!
# Match the roots and their complex conjugates
j = match_vectors(z1, Conj(z1))
print(data.frame(z = z1, j = j, `Conj(z[j])` = Conj(z1[j]), 
                 d = z1-Conj(z1[j])))
#>                          z j               Conj.z.j..
#> 1  0.1081845-3.805031e-23i 1  0.1081845+3.805031e-23i
#> 2  0.6907488+7.141664e-01i 3  0.6907488+7.141664e-01i
#> 3  0.6907488-7.141664e-01i 2  0.6907488-7.141664e-01i
#> 4 -1.9146323+7.502615e-15i 4 -1.9146323-7.502615e-15i
#> 5  7.3319114+1.630502e-15i 5  7.3319114-1.630502e-15i
#>                             d
#> 1  0.000000e+00-7.610062e-23i
#> 2  1.187939e-14-9.214851e-15i
#> 3 -1.187939e-14-9.214851e-15i
#> 4  0.000000e+00+1.500523e-14i
#> 5  0.000000e+00+3.261005e-15i
```
