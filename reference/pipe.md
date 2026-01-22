# Pipe operator

Re-export pipe operator `%>%` to turn function composition into a series
of imperative statements. For more extensive description, see function
`` `%>%` `` in package *magrittr*.

## Arguments

- lhs:

  First argument of the function of the right-hand-side of the pipe
  operator.

- rhs:

  Function whose first argument is given by the left-hand-side argument
  `lhs` of the pipe operator.

## Examples

``` r
x = array(stats::rnorm(2*1*3, sd = 0.01), dim = c(2,1,3))
# Instead of
pseries(polm(x))
#> ( 2 x 1 ) impulse response with maximum lag = 5 
#>        lag=0 [,1]    lag=1 [,1]   lag=2 [,1] lag=3 [,1] lag=4 [,1] lag=5 [,1]
#> [1,] -0.014836470  0.0204484548 -0.012872482          0          0          0
#> [2,] -0.006548665 -0.0001316217 -0.009101773          0          0          0
# you can write
x %>% polm() %>% pseries()
#> ( 2 x 1 ) impulse response with maximum lag = 5 
#>        lag=0 [,1]    lag=1 [,1]   lag=2 [,1] lag=3 [,1] lag=4 [,1] lag=5 [,1]
#> [1,] -0.014836470  0.0204484548 -0.012872482          0          0          0
#> [2,] -0.006548665 -0.0001316217 -0.009101773          0          0          0
```
