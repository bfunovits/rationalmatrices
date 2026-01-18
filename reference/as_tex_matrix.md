# TeX Matrix

TeX Matrix

## Usage

``` r
as_tex_matrix(x)
```

## Arguments

- x:

  matrix, where `paste(x[i,j])` returns a valid "TeX" string.

## Value

character string.

## Examples

``` r
as_tex_matrix(diag(1:2, nrow = 2, ncol = 3))
#> [1] "\\begin{pmatrix}\n  1 & 0 & 0  \\\\\n  0 & 2 & 0  \\\\\n\\end{pmatrix}"
```
