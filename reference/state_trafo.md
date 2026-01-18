# State Transformation

Applies a "state transformation" to a given rational matrix \\K(z) =
C(Iz^{-1} - A)^{-1}B +D\\ in statespace form. The parameter matrices are
transformed as \\A \rightarrow T A T^{-1}\\, \\B \rightarrow T B\\ and
\\C \rightarrow C T^{-1}\\, where \\T\\ is a non-singular (square)
matrix.

## Usage

``` r
state_trafo(obj, T, inverse = FALSE)
```

## Arguments

- obj:

  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object, represents a rational matrix \\K(z)\\.

- T:

  state transformation matrix. The parameter `T` must be a non-singular
  (s-by-s) matrix, or a vector of length s^2 (in this case T is coerced
  to an s-by-s matrix) or a vector of length s. In the latter case T is
  coerced to a diagonal s-by-s matrix.

- inverse:

  if TRUE, the transformation is reversed, i.e. \\A \rightarrow T^{-1} A
  T\\, \\B \rightarrow T^{-1} B\\ and \\C \rightarrow C T\\.

## Value

[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
object which represents the transformed state space realization.

## Examples

``` r
obj = stsp(A = c(0,0.2,1,-0.5),
           B = c(1,1), C = c(1,0))
           
# random state transformation
T = stats::rnorm(4)
obj1 = state_trafo(obj, T)

# obj and obj1 are equivalent, they produce the same IRF
testthat::expect_equivalent(pseries(obj), pseries(obj1))

# diagonal state transformation matrix
T = stats::rnorm(2)
obj1 = state_trafo(obj, T)

# revert the transformation
testthat::expect_equivalent(obj, state_trafo(obj1, T, inverse = TRUE))
if (FALSE) { # \dontrun{
state_trafo(obj, stats::rnorm(9)) # dimension of T does not fit
state_trafo(obj, c(1,0))          # T = diag(c(1,0)) is singular!
} # }
```
