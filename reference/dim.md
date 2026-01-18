# Dimensions of Objects

Retrieve the dimension and degrees of rational matrix objects. E.g. for
a polynomial matrix \\x(z)\\ (i.e. for a
[`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
object), `dim(x)` returns the vector `c(m,n,p)`, where `m,n` are the
respective number of rows and columns of the matrix and `p` is the
(polynomial) degree of the matrix. For a Laurent polynomial, there is an
additional element in the vector `c(m,n,p,min_deg)` returned by `dim(x)`
pertaining to the (possibly negative) minimal degree. For a rational
matrix in RMFD form \\x(z)=d(z)c^{-1}(z)\\ (i.e. for an
[`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
object), `dim(x)` returns the vector `c(m,n,p,q)`, where `m,n` are the
respective number of rows and columns of the matrix and `p,q` are the
(polynomial) degrees of the polynomial matrices \\c(z)\\ and \\d(z)\\
respectively.

## Usage

``` r
# S3 method for class 'polm'
dim(x)

# S3 method for class 'lpolm'
dim(x)

# S3 method for class 'lmfd'
dim(x)

# S3 method for class 'rmfd'
dim(x)

# S3 method for class 'stsp'
dim(x)

# S3 method for class 'pseries'
dim(x)

# S3 method for class 'zvalues'
dim(x)
```

## Arguments

- x:

  Object of type
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md),
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md),
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md),
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md),
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md),
  or
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md).

## Value

Returns a named vector of integers (`m,n` always refer to the number of
rows and columns of the rational matrix)

- `c(m,n,p,min_deg`:

  for Laurent polynomial matrices (i.e. for
  [`lpolm`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  objects).

- `c(m,n,p)`:

  for polynomial matrices (i.e. for
  [`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  objects).

- `c(m,n,p,q)`:

  for left matrix fraction descriptions (i.e. for
  [`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  objects).

- `c(m,n,p,q)`:

  for right matrix fraction descriptions (i.e. for
  [`rmfd`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
  objects).

- `c(m,n,s)`:

  for state space representations (i.e. for
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  objects).

- `c(m,n,lag.max)`:

  for power series expansions (i.e. for
  [`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  objects).

- `c(m,n,n.f)`:

  for frequency response functions (i.e. for
  [`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  objects).
