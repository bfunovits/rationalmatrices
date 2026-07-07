# Forward and Backward Bracket

For an
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object as input, `get_bwd` discards all coefficient matrices pertaining
to **negative** powers and returns a
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object with `min_deg = 0`. Similarly, `get_fwd` discards all coefficient
matrices pertaining to **non-negative** powers, and also returns an
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
object.

## Usage

``` r
get_fwd(lpolm_obj)

get_bwd(lpolm_obj)
```

## Arguments

- lpolm_obj:

  Laurent polynomial object
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)

## Value

Laurent polynomial object
[`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
without non-negative coefficients or without negative

## Details

Obtain the forward or backward part of a Laurent polynomial, i.e. apply
\\\[.\]\_-\\ or \\\[.\]\_+\\ to \$\$a(z) = a\_{-q} z^{-q} + \cdots +
a\_{-1} z^{-1} + a_0 + a_1 z^1 + \cdots + a_p z^p\$\$ and obtain for
`get_fwd` \$\$\[a(z)\]\_- = a\_{-q} z^{-q} + \cdots + a\_{-1} z^{-1}\$\$
or for `get_bwd` \$\$\[a(z)\]\_+ = a_0 + a_1 z^1 + \cdots + a_p z^p\$\$

## Examples

``` r
(lp = test_lpolm(degree_max = 2, degree_min = -2))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 2, and minimal degree >= -2
#>      z^-2 [,1] z^-1 [,1]  z^0 [,1]   z^1 [,1]  z^2 [,1]
#> [1,] -1.092631 -1.947076 -1.521579 -0.1541234 0.8862193
get_fwd(lp)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= -1, and minimal degree >= -2
#>      z^-2 [,1] z^-1 [,1]
#> [1,] -1.092631 -1.947076
(lp = test_lpolm(degree_max = 2, degree_min = -2))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 2, and minimal degree >= -2
#>      z^-2 [,1] z^-1 [,1]   z^0 [,1]  z^1 [,1] z^2 [,1]
#> [1,] -1.419259 0.3932093 -0.8518418 0.7281672 1.716373
get_bwd(lp)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 2, and minimal degree >= 0
#>        z^0 [,1]  z^1 [,1] z^2 [,1]
#> [1,] -0.8518418 0.7281672 1.716373

(lp = lpolm(1:3, min_deg = 2))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 4, and minimal degree >= 2
#>      z^2 [,1] z^3 [,1] z^4 [,1]
#> [1,]        1        2        3
get_bwd(lp)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 4, and minimal degree >= 0
#>      z^0 [,1] z^1 [,1] z^2 [,1] z^3 [,1] z^4 [,1]
#> [1,]        0        0        1        2        3

(lp = lpolm(1:3, min_deg = -1))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= -1
#>      z^-1 [,1] z^0 [,1] z^1 [,1]
#> [1,]         1        2        3
get_bwd(lp)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= 1, and minimal degree >= 0
#>      z^0 [,1] z^1 [,1]
#> [1,]        2        3

(lp = lpolm(1:3, min_deg = -5))
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= -3, and minimal degree >= -5
#>      z^-5 [,1] z^-4 [,1] z^-3 [,1]
#> [1,]         1         2         3
get_bwd(lp)
#> ( 1 x 1 ) Laurent polynomial matrix with degree <= -1, and minimal degree >= 0
```
