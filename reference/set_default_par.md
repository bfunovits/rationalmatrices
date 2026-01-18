# Set Default Graphical Parameters

Set default graphical parameters.

## Usage

``` r
set_default_par(m = 1, n = 1)
```

## Arguments

- m, n:

  integers. The procedure assumes that the device is split up into an
  m-by-n array of (sub-) figures and scales `cex.axis`, `cex.lab` and
  `cex.main` correspondingy.

## Value

invisible list with the original graphical parameters.

## Examples

``` r
opar = set_default_par(2,3)
print(opar)
#> $oma
#> [1] 0 0 0 0
#> 
#> $fig
#> [1] 0 1 0 1
#> 
#> $mar
#> [1] 5.1 4.1 4.1 2.1
#> 
#> $tcl
#> [1] -0.5
#> 
#> $mgp
#> [1] 3 1 0
#> 
#> $mex
#> [1] 1
#> 
#> $cex
#> [1] 1
#> 
#> $mex
#> [1] 1
#> 
#> $cex.axis
#> [1] 1
#> 
#> $cex.lab
#> [1] 1
#> 
#> $cex.main
#> [1] 1.2
#> 
#> $xaxs
#> [1] "r"
#> 
#> $yaxs
#> [1] "r"
#> 
#> $col.axis
#> [1] "black"
#> 
```
