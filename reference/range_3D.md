# Range of Values in a 3-D Array

`range_3D` computes the range of values of a three-dimensional array.

## Usage

``` r
range_3D(x, MARGIN = integer(0))
```

## Arguments

- x:

  `(m,n,k)`-dimensional numeric array

- MARGIN:

  integer(0), 1, 2 or (1:2). determines the "margin" over which to
  compute the range(s).

## Value

`(m,n,2)`-dimensional array.

## Details

If the input, `x` say, is an `(m,n,k)` dimensional array then `range_3D`
returns an `(m,n,2)` dimensional array, `m` say:

- `MARGIN = c(1,2)`: `m[i,j,]` contains the range of values in
  `x[i,j,]`.

- `MARGIN = 1`: `m[i,j,]` contains the range of values in `x[i,,]`.

- `MARGIN = 2`: `m[i,j,]` contains the range of values in `x[,j,]`.

- `MARGIN = integer(0)`: `m[i,j,]` contains the range of values in
  `x[,,]`.
