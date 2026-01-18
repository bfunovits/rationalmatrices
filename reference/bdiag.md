# Block Diagonal Matrix

Combine two or more matrices to a block-diagonal matrix. The functions
supports boolean, integer, numeric and complex matrices (and vectors).
The procedure makes some effort to retain the (col- / row-) names of the
inputs.

## Usage

``` r
bdiag(...)
```

## Arguments

- ...:

  matrices or vectors. Vectors are treated as diagonal matrices. If no
  input arguments are provided then `bdiag` returns `NULL`.

## Value

block diagonal matrix (or `NULL` if `bdiag` has been called without
inputs).

## Examples

``` r
A = matrix(TRUE, nrow = 2, ncol = 3)
colnames(A) = paste('A',1:3,sep ='.')
B = rep(2L,3)
names(B) = paste('B',1:3,sep = '.')
C = matrix(0, nrow = 3, ncol = 0)
rownames(C) = paste('C',1:3)
D = matrix(4, nrow = 2, ncol = 3)
E = matrix(complex(real=5), nrow = 2, ncol = 2)
rownames(E) = paste('E',1:2,sep='.')
colnames(E) = paste('E',1:2,sep='.')

bdiag()
#> NULL
bdiag(NULL,NULL)
#> NULL

X = bdiag(A,NULL) # NULL arguments are skipped
X
#>       A.1  A.2  A.3
#> [1,] TRUE TRUE TRUE
#> [2,] TRUE TRUE TRUE
str(X)            # output is of type 'logi'
#>  logi [1:2, 1:3] TRUE TRUE TRUE TRUE TRUE TRUE
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:3] "A.1" "A.2" "A.3"

X = bdiag(A,NULL,B)
X
#>     A.1 A.2 A.3 B.1 B.2 B.3
#>       1   1   1   0   0   0
#>       1   1   1   0   0   0
#> B.1   0   0   0   2   0   0
#> B.2   0   0   0   0   2   0
#> B.3   0   0   0   0   0   2
str(X)            # output is of type 'int'
#>  int [1:5, 1:6] 1 1 0 0 0 1 1 0 0 0 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:5] "" "" "B.1" "B.2" ...
#>   ..$ : chr [1:6] "A.1" "A.2" "A.3" "B.1" ...

# note the action of the "empty" (3 times 0) matrix C
X = bdiag(A,C,B)
X
#>     A.1 A.2 A.3 B.1 B.2 B.3
#>       1   1   1   0   0   0
#>       1   1   1   0   0   0
#> C 1   0   0   0   0   0   0
#> C 2   0   0   0   0   0   0
#> C 3   0   0   0   0   0   0
#> B.1   0   0   0   2   0   0
#> B.2   0   0   0   0   2   0
#> B.3   0   0   0   0   0   2
str(X)            # output is of type 'num'
#>  num [1:8, 1:6] 1 1 0 0 0 0 0 0 1 1 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:8] "" "" "C 1" "C 2" ...
#>   ..$ : chr [1:6] "A.1" "A.2" "A.3" "B.1" ...

X = bdiag(A,B,C,D,E)
X
#>      A.1  A.2  A.3  B.1  B.2  B.3                 E.1  E.2
#>     1+0i 1+0i 1+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#>     1+0i 1+0i 1+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#> B.1 0+0i 0+0i 0+0i 2+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#> B.2 0+0i 0+0i 0+0i 0+0i 2+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#> B.3 0+0i 0+0i 0+0i 0+0i 0+0i 2+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#> C 1 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#> C 2 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#> C 3 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i
#>     0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 4+0i 4+0i 4+0i 0+0i 0+0i
#>     0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 4+0i 4+0i 4+0i 0+0i 0+0i
#> E.1 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 5+0i 5+0i
#> E.2 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 0+0i 5+0i 5+0i
str(X)            # output is of type 'cplx'
#>  cplx [1:12, 1:11] 1+0i 1+0i 0+0i ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:12] "" "" "B.1" "B.2" ...
#>   ..$ : chr [1:11] "A.1" "A.2" "A.3" "B.1" ...

if (FALSE) { # \dontrun{
# the inputs must be vectors or matrices. arrays are not supported.
bdiag(A, array(1, dim = c(2,3,1)))

# character matrices are problematic, since it is not clear how to set the
# non diagonal elements. Therefore, the following statement throws an error.
bdiag(A, matrix('B', nrow = 2, ncol = 1), 3:4)
} # }
```
