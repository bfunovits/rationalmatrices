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
#> [1] "\\begin{pmatrix}\n  -1.4 - 1.4z & 0 & 0.4 + 0.4z  \\\\\n  1.4 + 1.4z & -0.1 - 0.1z & 0.1 + 0.1z  \\\\\n\\end{pmatrix}"
as_tex_matrixpoly(coefs, x = 'x', as_matrix_of_polynomials = FALSE)
#> [1] " \\begin{pmatrix}\n  -1.4 & 0 & 0.4  \\\\\n  1.4 & -0.1 & 0.1  \\\\\n\\end{pmatrix}  +\n \\begin{pmatrix}\n  -1.4 & 0 & 0.4  \\\\\n  1.4 & -0.1 & 0.1  \\\\\n\\end{pmatrix} x"
```
