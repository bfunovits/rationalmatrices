# Extract Parts of a Rational Matrix

The subsetting operation `x[,]` for rational matrices works analogously
to the subsetting of ordinary matrices. However, this operator is only
partly implemented for
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
and
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
objects. See the details below.  
The `$` operator may be used to extract the polynomial factors of a
left/right matrix fraction description
([`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
object). Furthermore one may retrieve the parameter matrices \\A,B,C,D\\
of a state space representation
([`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
object). For an
[`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
object, we may access the complex numbers at which the rational matrix
has been evaluated.

## Usage

``` r
# S3 method for class 'polm'
x[i, j]

# S3 method for class 'lpolm'
x[i, j]

# S3 method for class 'lmfd'
x[i, j]

# S3 method for class 'rmfd'
x[i, j]

# S3 method for class 'stsp'
x[i, j]

# S3 method for class 'pseries'
x[i, j]

# S3 method for class 'zvalues'
x[i, j]

# S3 method for class 'lmfd'
x$name

# S3 method for class 'rmfd'
x$name

# S3 method for class 'stsp'
x$name

# S3 method for class 'zvalues'
x$name
```

## Arguments

- x:

  a rational matrix, i.e. a
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  object.

- i, j:

  indices (integer or boolean vector)

- name:

  character: A,B,C,D for `stsp` objects, a,b for `lmfd` objects, c,d for
  `rmfd` objects and z,f for `zvalues` objects.

## Value

The subsetting operation `x[i,j]` returns a rational matrix of the same
class as the input `x`. The mode of the output of the `$` operator
depends on which "component" of the rational matrix is extracted. See
the details above.

## Details

- `x[]` or `x[,]` simply return the original object.

- `x[i]` returns a "vector", i.e. an \\(s,1)\\ dimensional matrix.

- `x[i,]`, `x[,j]` or `x[i,j]` return a rational matrix with rows
  selected by `i` and columns selected by `j`.

- **Note:** for
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  objects
  ([`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
  objects) only extraction of columns (respectively rows) is
  implemented. Therefore, e.g. `x[i,]` throws an error if `x` is an
  `lmfd` object.

- Note that "named" arguments are not supported (in order to simplify
  the coding).

- In order to have a finer control, one may e.g. use `unclass(x)[,,]`.

- `x$a`, `x$b` returns the left, respectively right factor of an LMFD
  \\x(z) = a^{-1}(z)b(z)\\. (`x` is an
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  object).

- `x$c`, `x$d` returns the right, respectively left factor of an RMFD
  \\x(z) = d(z)c^{-1}(z)\\. (`x` is an `rmfd` object.)

- `x$A`, `x$B`, `x$C`, `x$D` return the parameter matrices of the
  statespace realization \\x(z) = D + z (I- Az)^{-1}B\\. (`x` is an
  `stsp` object.)

- If `x` is an
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  object, then `x$z` returns the vector of complex points at which the
  rational matrix \\x\\ has been evaluated. Furthermore `x$f` gives the
  corresponding "frequencies", i.e. `x$f = -Arg(x$z)/(2*pi)`.

## Examples

``` r
# polynomial matrices 
a = test_polm(dim = c(3,2), degree = 1)
a[]          # returns a
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]      110   120      111   121
#> [2,]      210   220      211   221
#> [3,]      310   320      311   321
a[,]         # returns a
#> ( 3 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]      110   120      111   121
#> [2,]      210   220      211   221
#> [3,]      310   320      311   321
a[c(1,3,6)]  # returns a "vector" with the (1,1), (3,1) and (3,2) element of a
#> ( 3 x 1 ) matrix polynomial with degree <= 1 
#>      z^0 [,1] z^1 [,1]
#> [1,]      110      111
#> [2,]      310      311
#> [3,]      320      321
a[1,]        # returns the first row of a 
#> ( 1 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]      110   120      111   121
a[,2]        # returns the second column of a 
#> ( 3 x 1 ) matrix polynomial with degree <= 1 
#>      z^0 [,1] z^1 [,1]
#> [1,]      120      121
#> [2,]      220      221
#> [3,]      320      321
a[c(TRUE,FALSE,TRUE),c(FALSE, TRUE)] # returns a 2 by 1 matrix 
#> ( 2 x 1 ) matrix polynomial with degree <= 1 
#>      z^0 [,1] z^1 [,1]
#> [1,]      120      121
#> [2,]      320      321
a[c(1,1),c(2,1)] # returns a 2 by 2 matrix 
#> ( 2 x 2 ) matrix polynomial with degree <= 1 
#>      z^0 [,1]  [,2] z^1 [,1]  [,2]
#> [1,]      120   110      121   111
#> [2,]      120   110      121   111
# check with pseries/zvalues
all.equal(pseries(a[c(1,1),c(2,1)]), pseries(a)[c(1,1),c(2,1)])
#> [1] TRUE
all.equal(zvalues(a[c(1,1),c(2,1)]), zvalues(a)[c(1,1),c(2,1)])
#> [1] TRUE

if (FALSE) { # \dontrun{
a[i=1, j=2] # throws an error, since "named" arguments are not allowed.
} # }

# the subsetting operator [,] is only implemented for "lmfd" columns
(l = test_lmfd(dim = c(2,2), degrees = c(1,1)))
#> ( 2 x 2 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0 -0.2652126  0.1129953
#> [2,]        0     1 -1.0874728 -0.3687616
#> right factor b(z):
#>        z^0 [,1]       [,2]   z^1 [,1]       [,2]
#> [1,] -1.8138184 -0.3233435  0.5227212 -0.7382624
#> [2,] -0.7279944  0.1181781 -0.3947120  1.3197001
l[,1]
#> ( 2 x 1 ) left matrix fraction description a^(-1)(z) b(z) with degrees (p = 1, q = 1)
#> left factor a(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0 -0.2652126  0.1129953
#> [2,]        0     1 -1.0874728 -0.3687616
#> right factor b(z):
#>        z^0 [,1]   z^1 [,1]
#> [1,] -1.8138184  0.5227212
#> [2,] -0.7279944 -0.3947120

# the subsetting operator [,] is only implemented for "rmfd" rows
(r = test_rmfd(dim = c(2,2), degrees = c(1,1)))
#> ( 2 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>        z^0 [,1]       [,2]  z^1 [,1]        [,2]
#> [1,]  1.1587522  0.8627108 -1.636792  1.53937650
#> [2,] -0.6343578 -1.6815730  1.430290 -0.07721086
#> right factor c(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0 -0.5477175  0.2388897
#> [2,]        0     1  2.1165307 -0.5159694
r[1,]
#> ( 1 x 2 ) right matrix fraction description d(z) c^(-1)(z) with degrees deg(c(z)) = p = 1, deg(d(z)) = q = 1
#> left factor d(z):
#>      z^0 [,1]      [,2]  z^1 [,1]     [,2]
#> [1,] 1.158752 0.8627108 -1.636792 1.539376
#> right factor c(z):
#>      z^0 [,1]  [,2]   z^1 [,1]       [,2]
#> [1,]        1     0 -0.5477175  0.2388897
#> [2,]        0     1  2.1165307 -0.5159694
```
