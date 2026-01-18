# TeX Matrix Polynomial Filters

TeX Matrix Polynomial Filters

## Usage

``` r
as_tex_matrixfilter(coefs, x = "z", t = "t")
```

## Arguments

- coefs:

  3-dimensional array with the coefficients of the filter.

- x:

  (character string) polynomial variable or process variable.

- t:

  (character string) time/index variable.

## Value

character string.

## Examples

``` r
coefs = array(round(rnorm(2*3*1), 1), dim = c(2,3,2))

as_tex_matrixfilter(coefs, x = '\\epsilon', t = 's')
#> [1] " \\begin{pmatrix}\n  -0.1 & 0.1 & -0.5  \\\\\n  0.4 & -0.8 & -1.2  \\\\\n\\end{pmatrix} \\epsilon_{s} +\n \\begin{pmatrix}\n  -0.1 & 0.1 & -0.5  \\\\\n  0.4 & -0.8 & -1.2  \\\\\n\\end{pmatrix} \\epsilon_{s-1}"
```
