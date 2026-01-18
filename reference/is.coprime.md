# Left Prime and Left Coprime Polynomials

Check whether a polynomial is *left prime* or a pair of two polynomials
is *left coprime*. This check uses a (singular) pencil associated with
the polynomial(s). For more details see the vignette [Technical
Details](https://bfunovits.github.io/rationalmatrices/doc/technical_details.md).

## Usage

``` r
is.coprime(
  a,
  b = NULL,
  tol = sqrt(.Machine$double.eps),
  only.answer = TRUE,
  debug = FALSE
)
```

## Arguments

- a, b:

  If `a` is an
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  or an
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
  object, which represents a left matrix fraction description, say
  \\p^{-1}(z) q(z)\\, or a right MFD, say \\r(z) s^{-1}(z)\\, then the
  procedure tests whether the pair \\(p(z),q(z))\\ or
  \\(t(r(z)),t(q(z)))\\is left coprime.

  Otherwise the arguments `a` and `b` (if `b` is not NULL) must
  represent two compatible polynomial matrices, i.e. `a`, `b` must be
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects (or objects which may be coerced to `polm` objects). If `b` is
  NULL, the procedures checks whether \\a(z)\\ is left prime, otherwise
  the pair \\(a(z),b(z))\\ is checked for left coprimeness.

- tol:

  a tolerance parameter, which is used to decide the rank of certain
  matrices.

- only.answer:

  if TRUE, just return a logical (`TRUE` or `FALSE`). Otherwise a list
  with additional information is returned.

- debug:

  if TRUE, print some diagnostic information.

## Value

If `only.answer` is true then a logical (`TRUE` or `FALSE`) is returned.
Otherwise, a list with the following slots is returned. A more detailed
description of these items is given in the vignette [Technical
Details](https://bfunovits.github.io/rationalmatrices/doc/technical_details.md).

- answer:

  A boolean as above

- A,B:

  These matrices represent the pencil \\(A-Bz)\\ (in staircase form)
  which is used to check the left (co-)prime condition.

- m,n:

  Two integer vectors which code the structure of the staircase form.

- zeroes:

  If available, a vector of zeroes of the matrix \\(a(z),b(z))\\. If
  \\(a,b)\\ have no common zeroes (the left coprime case) then `zeroes`
  is an empty numeric vector. The case that \\(a(z),b(z))\\ is rank
  deficient for *all* \\z \in C\\ is coded with `z=NA`.

## Note

This procedure returns different objects, depending on the parameter
`only.answer`.

## Examples

``` r
# Ex 1: Two coprime polynomials ##################################################

# Generate two random (2 x 2) polynomial matrices with degree 2
set.seed(1803)
a = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)
b = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)

# Output: "only.answer = TRUE"
is.coprime(a, b, debug = FALSE, only.answer = TRUE)
#> [1] TRUE

# Output: "only.answer = FALSE"
out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
str(out)
#> List of 6
#>  $ answer: logi TRUE
#>  $ A     : num [1:6, 1:8] -1.67e-17 -1.00 0.00 0.00 0.00 ...
#>  $ B     : num [1:6, 1:8] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ zeroes: num(0) 
#>  $ m     : int [1:3] 2 2 2
#>  $ n     : num [1:3] 2 2 4
out$answer
#> [1] TRUE
out$zeroes
#> numeric(0)

# we could equivalently use the syntax: 
is.coprime(cbind(a,b))
#> [1] TRUE
is.coprime(lmfd(a,b))
#> [1] TRUE

# Ex 2: Two non-coprime polynomials with a finite number of common zeros #############
# Dimensions of a, b, and the common factor r
dim = 3
deg_aa = 1
deg_bb = 1
deg_r = 1

# Generate random polynomial matrices
a0 = a
b0 = b
# generate common factor 
r = test_polm(dim = c(2,2), degree = 1, random = TRUE, digits = 1)

# Generate polynomials with a common factor
a = r %r% a0
b = r %r% b0

out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
out$answer
#> [1] FALSE
out$zeroes
#> [1] -0.8401447  0.5105993


# Ex 3: Two non-coprime polynomials: Everywhere rank deficient ###################

# generate a common factor of rank 1 
r = test_polm(dim = c(2,1), degree = 1, random = TRUE, digits = 1) %r% 
    test_polm(dim = c(1,2), degree = 1, random = TRUE, digits = 1)

# Rank deficient matrices with common factor
a = r %r% a0
b = r %r% b0

out = is.coprime(a,b, only.answer = FALSE)
out$answer
#> [1] FALSE
out$zeroes
#> [1] NA


# Ex 4: Right-MFD ####

c = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)
d = test_polm(dim = c(2,2), degree = 2, random = TRUE, digits = 1)

# Output: "only.answer = TRUE"
is.coprime(t(c), t(d), debug = FALSE, only.answer = TRUE)
#> [1] TRUE

# Output: "only.answer = FALSE"
out = is.coprime(t(c), t(d), debug = FALSE, only.answer = FALSE)
str(out)
#> List of 6
#>  $ answer: logi TRUE
#>  $ A     : num [1:6, 1:8] -1.00 3.42e-17 0.00 0.00 0.00 ...
#>  $ B     : num [1:6, 1:8] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ zeroes: num(0) 
#>  $ m     : int [1:3] 2 2 2
#>  $ n     : num [1:3] 2 2 4
out$answer
#> [1] TRUE
out$zeroes
#> numeric(0)

# we could equivalently use the syntax: 
is.coprime(rbind(c,d))
#> [1] FALSE
is.coprime(rmfd(c,d))
#> [1] TRUE

# reset seed
set.seed(NULL)
```
