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
#> 1  0.1437411+9.484240e-01i 3  0.1437409+9.484245e-01i
#> 2 -1.1382540-1.859817e-18i 4 -1.1382546+8.997297e-17i
#> 3  0.1437411-9.484240e-01i 5  0.1437409-9.484245e-01i
#> 4  1.1040457-1.371466e+00i 2  1.1040450-1.371466e+00i
#> 5  1.1040457+1.371466e+00i 1  1.1040450+1.371466e+00i
#>                            d
#> 1 1.711940e-07-5.122472e-07i
#> 2 6.564261e-07-9.183279e-17i
#> 3 1.711940e-07+5.122472e-07i
#> 4 6.316591e-07-7.243216e-08i
#> 5 6.316591e-07+7.243216e-08i

# A polynomial with real coefficients has pairs of complex conjugate roots.
# However, the roots returned by "polyroot" in general do not have this 
# property!
# Match the roots and their complex conjugates
j = match_vectors(z1, Conj(z1))
print(data.frame(z = z1, j = j, `Conj(z[j])` = Conj(z1[j]), 
                 d = z1-Conj(z1[j])))
#>                          z j               Conj.z.j..
#> 1  0.1437411+9.484240e-01i 3  0.1437411+9.484240e-01i
#> 2 -1.1382540-1.859817e-18i 2 -1.1382540+1.859817e-18i
#> 3  0.1437411-9.484240e-01i 1  0.1437411-9.484240e-01i
#> 4  1.1040457-1.371466e+00i 5  1.1040457-1.371466e+00i
#> 5  1.1040457+1.371466e+00i 4  1.1040457+1.371466e+00i
#>                             d
#> 1 -2.220446e-16+0.000000e+00i
#> 2  0.000000e+00-3.719634e-18i
#> 3  2.220446e-16+0.000000e+00i
#> 4 -8.881784e-16+0.000000e+00i
#> 5  8.881784e-16+0.000000e+00i
```
