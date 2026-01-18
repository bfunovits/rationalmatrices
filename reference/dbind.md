# Bind Arrays

`dbind(d, x, y, ...)` concatenates/binds an arbitrary number of arrays
`x, y,...` along the dimension `d`. For matrices,
`dbind(d = 1, x, y, ...)` is (essentially) equivalent to
[`rbind`](https://rdrr.io/r/base/cbind.html) and
`dbind(d = 2, x, y, ...)` corresponds to
[`cbind`](https://rdrr.io/r/base/cbind.html). If the number of
dimensions of an argument, `x` say, is less than `d`, then this argument
is treated as an array of dimension `c(dim(x),1,..,1)`.

## Usage

``` r
dbind(d = 1, ...)
```

## Arguments

- d:

  Integer. Concatenate arrays along the dimension `d`.

- ...:

  arrays.

## Value

Array

## Details

The procedure makes some effort to keep the `dimnames` attribute of the
arguments.

## Examples

``` r
x = test_array(dim = c(2,3,1))
y = test_array(dim = c(2,3,1), dimnames = TRUE)
z = test_array(dim = c(2,3,3), dimnames = TRUE)

# Bind along dimension 1 (row-binding for matrices)
dbind(1, x)
#> , , 1
#> 
#>      [,1] [,2] [,3]
#> [1,]  111  121  131
#> [2,]  211  221  231
#> 
dbind(1, x, y)
#> , , C = C=1
#> 
#>       B
#> A      B=1 B=2 B=3
#>   [1,] 111 121 131
#>   [2,] 211 221 231
#>   [3,] 111 121 131
#>   [4,] 211 221 231
#> 

# Bind along dimension 2 (col-binding for matrices)
dbind(2, x, y)
#> , , C = C=1
#> 
#>      B
#> A     [,1] [,2] [,3] [,4] [,5] [,6]
#>   A=1  111  121  131  111  121  131
#>   A=2  211  221  231  211  221  231
#> 

# Bind along dimension 3
dbind(3, x, y)
#> , , 1
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 111 121 131
#>   A=2 211 221 231
#> 
#> , , 2
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 111 121 131
#>   A=2 211 221 231
#> 
dbind(3, x, y, z)
#> , , 1
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 111 121 131
#>   A=2 211 221 231
#> 
#> , , 2
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 111 121 131
#>   A=2 211 221 231
#> 
#> , , 3
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 111 121 131
#>   A=2 211 221 231
#> 
#> , , 4
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 112 122 132
#>   A=2 212 222 232
#> 
#> , , 5
#> 
#>      B
#> A     B=1 B=2 B=3
#>   A=1 113 123 133
#>   A=2 213 223 233
#> 

# Example that throws an error
if (FALSE) { # \dontrun{
dbind(1, x, y, z) # throws an error, since the array x,y,z are not compatible
} # }
```
