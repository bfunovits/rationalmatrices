# Check the Minimality of a Statespace Realization

Check whether a statespace realization of a rational matrix is minimal.
The procedure constructs the Hankel matrix of the impulse response
coefficients with \\s\\ block rows and \\s\\ block columns, where \\s\\
is the statespace dimension of the given statespace realization. If this
statespace realization is minimal then the Hankel matrix has rank \\s\\.
Therefore the procedure returns `TRUE` if the \\s\\-th singular values
of the Hankel matrix is larger than `tol`.

## Usage

``` r
is.minimal(x, ...)

# S3 method for class 'stsp'
is.minimal(x, tol = sqrt(.Machine$double.eps), only.answer = TRUE, ...)
```

## Arguments

- x:

  ([`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object) rational matrix in statespace form.

- ...:

  not used.

- tol:

  a tolerance parameter, which is used to decide the rank of the Hankel
  matrix of the impulse response coefficients.

- only.answer:

  if TRUE, just return a logical (TRUE or FALSE). Otherwise a list with
  additional information is returned.

## Value

If `only.answer` is true then a logical (`TRUE` or `FALSE`) is returned.
Otherwise, a list with the following slots is returned.

- answer:

  A boolean as above.

- H:

  Hankel matrix of impulse response coefficients (with s block rows and
  s block columns).

- sv:

  The singular values of H.

- s0:

  (integer) (estimate of the) rank of H, i.e. the minimal statespace
  dimension.

## Details

The procedure does not check whether the statespace realization is
*observable* and/or *controllable*. To this end one may compute the
observability/controllability matrices
([`obs_matrix`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md),
[`ctr_matrix`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md))
or the corresponding Grammians
([`grammians`](https://bfunovits.github.io/rationalmatrices/reference/grammians.md)).

## Note

This procedure returns different objects, depending on the parameter
`only.answer`.

## See also

[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
and
[`bhankel`](https://bfunovits.github.io/rationalmatrices/reference/bhankel.md).

## Examples

``` r
x = test_stsp(dim = c(2,2), s = 2)
is.minimal(x)
#> [1] TRUE
# note that operations on "stsp" objects may return non minimal realizations
# 
is.minimal(rbind(x, x), only.answer = FALSE)[c('answer','sv','s0')]
#> $answer
#> [1] FALSE
#> 
#> $sv
#> [1] 1.688895e+01 1.399243e+00 1.405004e-15 3.567018e-16 2.183431e-16
#> [6] 2.817023e-17 2.006397e-17 5.892006e-18
#> 
#> $s0
#> [1] 2
#> 
is.minimal(x %r% (x^(-1)), only.answer = FALSE)[c('answer','sv','s0')]
#> $answer
#> [1] FALSE
#> 
#> $sv
#> [1] 9.143145e-15 2.289160e-15 1.002634e-15 4.436468e-16 2.877783e-16
#> [6] 2.066874e-16 6.533035e-17 2.352831e-17
#> 
#> $s0
#> [1] 0
#> 

is.minimal(test_stsp(dim = c(2,0), s = 2))
#> [1] FALSE
is.minimal(test_stsp(dim = c(0,2), s = 0))
#> [1] TRUE
```
