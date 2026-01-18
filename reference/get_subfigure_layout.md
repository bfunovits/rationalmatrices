# Subfigure Layout

The device is split up into an (m-by-n) array of (sub-) figures. This
tool computes the corresponding figure region and margins for each of
the (sub-) figures.

## Usage

``` r
get_subfigure_layout(margins)
```

## Arguments

- margins:

  A named list with slots `bottom`, `left`, `top` and `right` which
  determine the margins of the sub-figures. E.g. `bottom` is an
  `m`-dimensional vector where `bottom[i]` gives the bottom margin (in
  "line units") of the sub-figures in the `i`-th row.

## Value

A list with slots

- omi:

  4-dimensional vector with the "outer margin" in inches.

- mar:

  (m,n,4)-dimensional array, where `mar[i,j,]` contains the margins (in
  "line" units) of the (i,j)-th sub-figure.

- fig:

  (m,n,4)-dimensional array, where `fig[i,j,]` contains the coordinates
  (in "NDC" units) of the figure region of the (i,j)-th sub-figure.

## See also

[`par`](https://rdrr.io/r/graphics/par.html) for more details.
