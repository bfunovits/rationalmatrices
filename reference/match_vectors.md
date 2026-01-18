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
#> 1  0.7876635+0.000000e+00i 2  0.7876625-6.462349e-27i
#> 2 -0.7046533-4.061788e-20i 1 -0.7046530-4.063111e-20i
#> 3  0.2719756-1.254518e+00i 5  0.2719748-1.254519e+00i
#> 4  0.2719756+1.254518e+00i 4  0.2719748+1.254519e+00i
#> 5  2.4999558-3.323473e-16i 3  2.4999524-3.323469e-16i
#>                             d
#> 1  9.780530e-07+6.462349e-27i
#> 2 -3.806902e-07+1.323489e-23i
#> 3  8.392146e-07+5.835721e-07i
#> 4  8.392146e-07-5.835721e-07i
#> 5  3.420185e-06-3.363525e-22i

# A polynomial with real coefficients has pairs of complex conjugate roots.
# However, the roots returned by "polyroot" in general do not have this 
# property!
# Match the roots and their complex conjugates
j = match_vectors(z1, Conj(z1))
print(data.frame(z = z1, j = j, `Conj(z[j])` = Conj(z1[j]), 
                 d = z1-Conj(z1[j])))
#>                          z j               Conj.z.j..
#> 1  0.7876635+0.000000e+00i 1  0.7876635+0.000000e+00i
#> 2 -0.7046533-4.061788e-20i 2 -0.7046533+4.061788e-20i
#> 3  0.2719756-1.254518e+00i 4  0.2719756-1.254518e+00i
#> 4  0.2719756+1.254518e+00i 3  0.2719756+1.254518e+00i
#> 5  2.4999558-3.323473e-16i 5  2.4999558+3.323473e-16i
#>                             d
#> 1  0.000000e+00+0.000000e+00i
#> 2  0.000000e+00-8.123575e-20i
#> 3  5.551115e-17+2.220446e-16i
#> 4 -5.551115e-17+2.220446e-16i
#> 5  0.000000e+00-6.646946e-16i
```
