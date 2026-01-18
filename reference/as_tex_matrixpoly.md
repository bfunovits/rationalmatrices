# TeX Matrix Polynomials

TeX Matrix Polynomials

## Usage

``` r
as_tex_matrixpoly(coefs, x = "z", as_matrix_of_polynomials = TRUE)
```

## Arguments

- coefs:

  3-dimensional array with the coefficients of the matrix polynomial.

- x:

  (character string) polynomial variable or process variable.

- as_matrix_of_polynomials:

  boolean.

## Value

character string.

## Examples

``` r
coefs = array(round(rnorm(2*3*1), 1), dim = c(2,3,2))

as_tex_matrixpoly(coefs)
#> [1] "\\begin{pmatrix}\n  1.5 + 1.5z & 0.3 + 0.3z & -0.6 - 0.6z  \\\\\n  -0.8 - 0.8z & -0.5 - 0.5z & 1.5 + 1.5z  \\\\\n\\end{pmatrix}"
as_tex_matrixpoly(coefs, x = 'x', as_matrix_of_polynomials = FALSE)
#> [1] " \\begin{pmatrix}\n  1.5 & 0.3 & -0.6  \\\\\n  -0.8 & -0.5 & 1.5  \\\\\n\\end{pmatrix}  +\n \\begin{pmatrix}\n  1.5 & 0.3 & -0.6  \\\\\n  -0.8 & -0.5 & 1.5  \\\\\n\\end{pmatrix} x"
```
