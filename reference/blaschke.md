# Blaschke Factors

The Blaschke factor at \\\alpha\\ is the rational function \$\$B(z) :=
\frac{1-\bar{\alpha}z}{-\alpha + z}\$\$ This is an all-pass function
with a pole at \\z=\alpha\\ and a zero at \\z=1/\bar{\alpha}\\. The
function `blaschke(alpha)` returns this rational \\(1 \times 1)\\ matrix
in
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
form. Clearly \\B(z)\\ has complex coefficients, if \\\alpha\\ is
complex.  
The call `blaschke2(alpha, row=NULL)` computes the product of the
Blaschke factors at \\\alpha\\ and at \\\bar{\alpha}\\, i.e. the
rational function \$\$B\_{s}(z) := \frac{1-2\Re(\alpha)z + \|\alpha\|^2
z^2}{\|\alpha\|^2 -2\Re(\alpha) z + z^2}\$\$  
If `blaschke2` is called with an optional argument `w` (a non zero
complex vector of length 2) then `blaschke2` constructs a \\(2 \times
2)\\ rational, all-pass matrix of the form \$\$B\_{2}(z) := a^{-1}(z)
b(z)\$\$ where \\a(z), b(z)\\ are two \\(2 \times 2)\\ polynomial
matrices (with real coefficients) of degree one. This matrix is
constructed such that the column space of \\a(\alpha)\\ is spanned by
\\\bar{w}\\ and the column space of \\a(\bar{\alpha})\\ is spanned by
the vector \\w\\.

## Usage

``` r
blaschke(alpha)

blaschke2(alpha, w = NULL, tol = 100 * .Machine$double.eps)
```

## Arguments

- alpha:

  complex or real scalar, represents \\\alpha\\.

- w:

  `NULL` or a (complex) vector of length 2.

- tol:

  Tolerance (used to decide whether `alpha` has modulus equal to one).

## Value

[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
object, which represents the constructed "Blaschke factors" \\B(z)\\,
\\B_s(z)\\ or \\B_2(z)\\.

## Note

The routine `blaschke2` throws an error if \\\alpha\\ is not complex
(i.e. the imaginary part is zero). If \\\alpha\\ is close to the unit
circle then `blaschke2(alpha, w)` simply returns an
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
representation of the bivariate identity matrix. If \\w\\ and
\\\bar{w}\\ are almost linearly dependent, then an error is thrown.

## Examples

``` r
# Blaschke factor with a real alpha
(B = blaschke(1.5))
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1] z^1 [,1]
#> [1,]     -1.5        1
#> right factor b(z):
#>      z^0 [,1] z^1 [,1]
#> [1,]        1     -1.5
zvalues(B) %>% abs()
#> ( 1 x 1 ) frequency response
#>      z[1] [,1] z[2] [,1] z[3] [,1] z[4] [,1] z[5] [,1]
#> [1,]         1         1         1         1         1

# Blaschke factor with a complex alpha
(B = blaschke(complex(real = 1.5, imaginary = 0.5)))
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>       z^0 [,1] z^1 [,1]
#> [1,] -1.5-0.5i     1+0i
#> right factor b(z):
#>      z^0 [,1]  z^1 [,1]
#> [1,]     1+0i -1.5+0.5i
zvalues(B) %>% abs()
#> ( 1 x 1 ) frequency response
#>      z[1] [,1] z[2] [,1] z[3] [,1] z[4] [,1] z[5] [,1]
#> [1,]         1         1         1         1         1

# product of the Blaschke factors at alpha and Conj(alpha) 
# this gives a scalar, rational, all-pass matrix with real coefficients
(B = blaschke2(complex(real = 1.5, imaginary = 0.5)))
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 2)
#> left factor a(z):
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]      2.5       -3        1
#> right factor b(z):
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]        1       -3      2.5
zvalues(B) %>% abs()
#> ( 1 x 1 ) frequency response
#>      z[1] [,1] z[2] [,1] z[3] [,1] z[4] [,1] z[5] [,1]
#> [1,]         1         1         1         1         1

#############################################################
# a "bivariate" Blaschke factor 

# case 1: alpha is "outside the unit circle" ################
(alpha = complex(real = 1.5, imaginary = 0.5))
#> [1] 1.5+0.5i
(w = complex(real = c(0.1,0.9), imaginary = c(0.75,-0.5)))
#> [1] 0.1+0.75i 0.9-0.50i

(B = blaschke2(alpha, w = w))
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0 -0.6786207 -0.1579310
#> [2,]        0     1  0.2924138 -0.5213793
#> right factor b(z):
#>       z^0 [,1]       [,2]  z^1 [,1]      [,2]
#> [1,] 0.6616180 -0.2854763 -0.976275  0.114544
#> [2,] 0.2214899  0.5090096  0.000000 -1.024302
# B(z) is all-pass 
print(zvalues(B) %r% Ht(zvalues(B)), digits = 3)
#> ( 2 x 2 ) frequency response
#>      z[1] [,1]  [,2] z[2] [,1]  [,2] z[3] [,1]  [,2] z[4] [,1]  [,2] z[5] [,1]
#> [1,]      1+0i  0+0i      1+0i  0+0i      1+0i  0+0i      1+0i  0+0i      1+0i
#> [2,]      0+0i  1+0i      0+0i  1+0i      0+0i  1+0i      0+0i  1+0i      0+0i
#>       [,2]
#> [1,]  0+0i
#> [2,]  1+0i

# B(z) has poles at z=alpha, z=Conj(alpha) and 
# zeroes at z=1/alpha and z=1/Conj(alpha)
poles(B)
#> [1] 1.5-0.5i 1.5+0.5i
zeroes(B)
#> [1] 0.6-0.2i 0.6+0.2i

# The column space of a(alpha) is spanned by the vector Conj(w).
max(abs( Conj(c(-w[2], w[1])) %*% zvalue(B$a, alpha) ))
#> [1] 1.475229e-16

# case 2: alpha is "inside the unit circle" #################
(alpha = 1 / complex(real = 1.5, imaginary = 0.5))
#> [1] 0.6-0.2i
(w = complex(real = c(0.1,0.9), imaginary = c(0.75,-0.5)))
#> [1] 0.1+0.75i 0.9-0.50i

(B = blaschke2(alpha, w = w))
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>        z^0 [,1]       [,2] z^1 [,1]  [,2]
#> [1,] -0.6786207 -0.1579310        1     0
#> [2,]  0.2924138 -0.5213793        0     1
#> right factor b(z):
#>      z^0 [,1]      [,2]   z^1 [,1]       [,2]
#> [1,] 0.976275 -0.114544 -0.6616180  0.2854763
#> [2,] 0.000000  1.024302 -0.2214899 -0.5090096
# B(z) is all-pass 
print(zvalues(B) %r% Ht(zvalues(B)), digits = 3)
#> ( 2 x 2 ) frequency response
#>      z[1] [,1]  [,2] z[2] [,1]  [,2] z[3] [,1]  [,2] z[4] [,1]  [,2] z[5] [,1]
#> [1,]      1+0i  0+0i      1+0i  0+0i      1+0i  0+0i      1+0i  0+0i      1+0i
#> [2,]      0+0i  1+0i      0+0i  1+0i      0+0i  1+0i      0+0i  1+0i      0+0i
#>       [,2]
#> [1,]  0+0i
#> [2,]  1+0i

# B(z) has poles at z=alpha, z=Conj(alpha) and 
# zeroes at z=1/alpha and z=1/Conj(alpha)
poles(B)
#> [1] 0.6-0.2i 0.6+0.2i
zeroes(B)
#> [1] 1.5-0.5i 1.5+0.5i

# The column space of a(alpha) is spanned by the vector Conj(w).
max(abs( Conj(c(-w[2], w[1])) %*% zvalue(B$a, alpha) ))
#> [1] 1.494683e-16

# case 3: alpha is "on the unit circle" #####################
alpha = alpha / Mod(alpha)
blaschke2(alpha) %>% print(digits = 2)
#> ( 1 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 2, q = 2)
#> left factor a(z):
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]        1     -1.9        1
#> right factor b(z):
#>      z^0 [,1] z^1 [,1] z^2 [,1]
#> [1,]        1     -1.9        1
blaschke2(alpha, w = w) %>% print(digits = 2)
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 0, q = 0)
#> left factor a(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
#> right factor b(z):
#>      z^0 [,1]  [,2]
#> [1,]        1     0
#> [2,]        0     1
```
