# Sequence generation

As opposed to the standard [`seq`](https://rdrr.io/r/base/seq.html)
command the function `iseq` returns an empty vector if the starting
value is larger than the end value.

## Usage

``` r
iseq(from = 1, to = 1)
```

## Arguments

- from, to:

  the starting and end values of the sequence.

## Value

`iseq` returns an empty integer vector if `to` is less than `from` and
`seq(from,to)` else.

## Details

More general than [seq_len](https://rdrr.io/r/base/seq.html) because the
sequence does not need to start from 1.

## Examples

``` r
iseq(0,1) # => c(0,1)
#> [1] 0 1
iseq(1,0) # => integer(0)
#> integer(0)
seq(1,0)  # => c(1,0)
#> [1] 1 0
```
