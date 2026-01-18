# Style Parameters

The plotting functions use a list of common style parameters.

## Usage

``` r
style_parameters(style = c("bw", "bw2", "gray", "colored"))
```

## Arguments

- style:

  character indicating the desired "style".

## Value

A list with slots:

- border.grid:

  background color of the plot region

- col.grid,lty.grid,lwd.grid:

  color, line style and line width of the grid lines

- col.axis:

  axis color

- col.box,lwd.box:

  color and linew idth of the box surrounding the plot region

- bg.labels,border.labels:

  background and border color for axis labels

- col.labels:

  text color of the axis labels

## Examples

``` r
style_parameters('colored')
#> $bg.grid
#> [1] "#E6E6E6"
#> 
#> $border.grid
#> [1] NA
#> 
#> $col.grid
#> [1] "white"
#> 
#> $lty.grid
#> [1] "solid"
#> 
#> $lwd.grid
#> [1] 1
#> 
#> $col.axis
#> [1] "#808080"
#> 
#> $col.box
#> [1] NA
#> 
#> $lwd.box
#> [1] 1
#> 
#> $bg.labels
#> [1] "orange"
#> 
#> $border.labels
#> [1] NA
#> 
#> $col.labels
#> [1] "black"
#> 
```
