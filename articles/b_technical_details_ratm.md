# Rational Matrices: Technical Details

## Introduction

This vignette provides technical implementation details and mathematical
theory underlying the `rationalmatrices` package. It complements the
main user guide (Rational Matrices vignette) with in-depth coverage of
algorithms, numerical methods, and theoretical foundations.

**Topics covered:** - **Left Coprime Polynomials** - Part of Analysis &
Properties (#5) - **Column Reduced Matrix** - Supporting theory for
Normal Forms (#4) - **Realization Algorithms** - Deep dive into
Realization Algorithms (#2) - **Operations with State-Space
Realizations** - Theory for State-Space Tools (#6) - **Schur, QZ,
Lyapunov** - Numerical Methods (#7) - **Reflect Poles and Zeroes** -
Theory for Pole/Zero Reflection (#8)

For implementation details of these topics, see the source files listed
in CLAUDE.md and the main vignette.

------------------------------------------------------------------------

## Left Coprime Polynomials

### Functional Area: Analysis & Properties (#5)

See `R/is_methods.R` for implementation.

A polynomial matrix $c$ is called *left prime*, if $c(z)$ has full row
rank everywhere in the complex plane. Clearly this implies that $c$ is
square or “wide”, i.e. if $c$ is $(m \times n)$-dimensional then
$m \leq n$ must hold.

A pair $(a,b)$ of (compatible) polynomial matrices is called *left
coprime* if the matrix $c = \lbrack a,b\rbrack$ is left prime. This case
is important for the structure of left matrix fraction descriptions.
Suppose $c(z) = a^{- 1}b(z)$, where $a$ is a square, non singular
polynomial matrix. If the pair is $(a,b)$ is *not* left coprime, then we
may cancel a common, non unimodular, factor and thus obtain a “simpler”
representation for $c(z)$.

We discuss two strategies for testing whether a polynomial matrix is
left prime. The first approach uses a Hermite normal form (HNF) and the
second one tackles this problem via a (singular) pencil.

The approach via singular pencils is numerically stable and thus is
implemented in the function
[`is.coprime()`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md).
See also the examples below.

First note that the problem is easy to tackle in the case that $c$ has
degree zero. Here we just have to check the left kernel of the
coefficient matrix of $z^{0}$. Also the case where $c$ is “tall”,
i.e. the number of rows is larger than the number of columns is trivial.
In this case $c$ clearly cannot be left prime, since the rank of $c(z)$
is smaller than the number of rows for any $z \in {\mathbb{C}}$.

Therefore throughout this section, we assume that $c$ is an
$(m \times n)$-dimensional polynomial matrix with degree $p > 0$ and
$m \leq n$.

### Hermite Normal Form

Consider the *row Hermite form* obtained from elementary column
operations $$c(z) = h(z)u(z) = h_{1}(z)u_{1}(z)$$ where $u$ is an
$(n \times n)$-dimensional unimodular modular matrix, $h$ is an
$(m \times n)$-dimensional lower (quasi) triangular matrix.

The $(m \times m)$-dimensional matrix $h_{1}$ is obtained from $h$ by
selecting the first $m$ columns and $u_{1}$ consists of the first $m$
rows of $u$. Note that $h_{1}$ is a square, lower triangular matrix,
where the diagonal elements are either monic or zero, and if a diagonal
element is non zero, then all elements to the left of this element have
a lower degree.

Furthermore note that the rank of $c(z)$ is equal to the rank of
$h_{1}(z)$ for all $z \in {\mathbb{C}}$. Therefore the rank of $c(z)$ is
less than $m$ if and only $z$ is a zero of the product of the diagonal
entries: $h_{11}h_{22}\cdots h_{mm}$. The following cases may occur

1.  All diagonal elements $h_{ii}$ are constant (and hence
    $h_{1} = I_{m}$ must hold). In this case
    $\text{rk}\left( c(z) \right) = m$ for all $z \in {\mathbb{C}}$ and
    thus $c$ is *coprime*.
2.  One of the diagonal elements is zero: In this case
    $\text{rk}\left( c(z) \right) < m$ holds for all
    $z \in {\mathbb{C}}$ and thus $c$ is *not coprime*.
3.  The diagonal elements are non zero and at least one of them has a
    degree larger than zero. In this case
    $\text{rk}\left( c(z) \right) = m$ except for the zeroes of the
    polynomial $\left( h_{11}h_{22}\cdots h_{mm} \right)$. Also in this
    case, $c$ is *not coprime*.

Hence, the polynomial $c$ is left prime if and only if $h_{1} = I_{m}$.

### Generalized Pencils

As reference for generalized pencils, see

- \[@Gantmacher1959, Chapter XII (in particular §4 on page 35ff)\]
- \[@Kailath80, page 393ff\]
- \[@DemmelKagstroem1993\]

We first construct a *pencil* $P(z) = (A - Bz)$ such that for any
$z \in {\mathbb{C}}$ the left kernel of $c(z)$ and of $P(z)$ have the
same dimension. This means that $c$ is left prime if and only if $P$ is
left prime. Furthermore the zeroes of $P$ and of $c$ are the same.

The condition
$$uc(z) = u\left( c_{0} + c_{1}z + \cdots + c_{p}z^{p} \right) = 0$$ for
$z \in {\mathbb{C}}$ and $u \in {\mathbb{C}}^{1 \times m},u \neq 0$ is
equivalent to the condition that
$$\left( v_{1},v_{2},\ldots,v_{p} \right)\begin{pmatrix}
c_{0} & 0 & \cdots & 0 \\
0 & I_{n} & \cdots & 0 \\
0 & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & I_{n} \\
 & & & 
\end{pmatrix} = \left( v_{1},v_{2},\ldots,v_{p} \right)\begin{pmatrix}
{- c_{1}} & \cdots & {- c_{p - 1}} & {- c_{p}} \\
I_{n} & \cdots & 0 & 0 \\
\vdots & \ddots & \vdots & \vdots \\
0 & \cdots & I_{n} & 0
\end{pmatrix}z$$ for $v_{1} \in {\mathbb{C}}^{1 \times m}$,
$v_{1} \neq 0$ and $v_{j} \in {\mathbb{C}}^{1 \times n}$ for
$i = 2,\ldots,p$. This follows immediately be rewriting the above
condition as $$\begin{array}{rcl}
{v_{2}z} & = & {v_{1}c_{0} + v_{1}c_{1}z = v_{1}\left( c_{0} + c_{1}z \right)} \\
{v_{3}z^{2}} & = & {\left( v_{2} + v_{1}c_{2}z \right)z = v_{1}\left( c_{0} + c_{1}z + c_{2}z^{2} \right)} \\
{v_{4}z^{3}} & = & {\left( v_{3} + v_{1}c_{3}z \right)z^{2} = v_{1}\left( c_{0} + c_{1}z + c_{2}z^{2} + c_{3}z^{3} \right)} \\
 & \vdots & \\
{v_{p}z^{p - 1}} & = & {v_{1}\left( c_{0} + c_{1}z + \cdots + c_{p - 1}z^{p - 1} \right)} \\
{- v_{1}c_{p}z^{p}} & = & {v_{p}z^{p - 1}}
\end{array}$$ see also @Kailath80, page 393ff. Combining the last two
equations gives
$$0 = v_{p}z^{p - 1} + v_{1}c_{p}z^{p} = v_{1}\left( c_{0} + c_{1}z + \cdots + c_{p}z^{p} \right).$$

Therefore the matrix $c(z)$ has a non trivial left kernel if and only if
the pencil $P(z) = (A - Bz)$,  
$$A = \begin{pmatrix}
c_{0} & 0 & \cdots & 0 \\
0 & I_{n} & \cdots & 0 \\
0 & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & I_{n} \\
 & & & 
\end{pmatrix},\quad B = \begin{pmatrix}
{- c_{1}} & \cdots & {- c_{p - 1}} & {- c_{p}} \\
I_{n} & \cdots & 0 & 0 \\
\vdots & \ddots & \vdots & \vdots \\
0 & \cdots & I_{n} & 0
\end{pmatrix},$$ has a non trivial left kernel. This implies $c$ is left
prime if and only if $P$ is left prime.

#### Staircase Form

In order to analyze the left kernel this pencil, we transform the pencil
(by multiplication with orthogonal matrices from the left and the right)
into a so called *staircase form*. To simplify notation we set
$M = \left( m + n(p - 1) \right)$ and $N = np$, i.e. $(M,N)$ now denotes
the dimension of the pencil $P(z)$. Note that $M - N = m - n \neq 0$ in
general, which means that we have to deal with *non regular* pencils.

The procedure described below works for an arbitary pencil, i.e. we do
not impose the above special structure and we also consider the case
$M > N$.

1.  If the pencil $P(z) = (A - Bz)$ is “tall”, i.e. $M > N$, then $P(z)$
    has a non trivial left kernel for all $z \in {\mathbb{C}}$.  
    In this case $P$ (and the matrix $c$) is *not* left prime.

2.  If the pencil $P(z) = (A - Bz)$ is square, i.e. $M = N$, and $B$ is
    non singular, then the zeroes of the pencil $P(z)$ (and the zeroes
    of $c(z)$) are the eigenvalues of the matrix $AB^{- 1}$.  
    In this case $P$ (and the matrix $c$) is *not* left prime.

Now suppose that the *right kernel* of $B$ is $N_{1} > 0$ dimensional.
If $N_{1} = N$, i.e. if $B = 0$, then two cases may occur:

3.  The matrix $A$ has full row rank.  
    In this case the pencil $P$ (and the matrix $c$) is left prime.
4.  The matrix $A$ does not have full row rank and hence $P(z)$ has a
    non trivial left kernel for all $z \in {\mathbb{C}}$.  
    In this case $P$ (and the matrix $c$) is *not* left prime.

Finally suppose that $0 < N_{1} < N$ holds and let
$0 \leq M_{1} \leq N_{1}$ denote the rank of the first $N_{1}$ columns
of $AV$. Then there exists an orthogonal transformation, $U$ say, such
that the first $M_{1}$ rows of the first $N_{1}$ columns of $U\prime AV$
have full rank $M_{1}$ and the remaining rows are zero.

Using the same symbols for the transformed pencil we get
$$P(z) = \begin{pmatrix}
A_{11} & A_{12} \\
0_{{(M - M_{1})} \times N_{1}} & A_{22}
\end{pmatrix} - \begin{pmatrix}
0_{M_{1} \times N_{1}} & B_{12} \\
0_{{(M - M_{1})} \times N_{1}} & B_{22}
\end{pmatrix}z,\quad A_{11} \in {\mathbb{R}}^{M_{1} \times N_{1}}.$$
Clearly $P(z)$ has full row rank if and only if the *reduced* pencil
$\left( A_{22} - B_{22}z \right)$ has full row rank.

Now we iterate this procedure with the reduced pencil
$\left( A_{22} - B_{22}z \right)$ until we end up with the following
staircase form of the pencil: $$P(z) = \begin{pmatrix}
A_{11} & \cdots & A_{1,k - 1} & A_{1k} \\
\vdots & \ddots & \vdots & \vdots \\
0_{M_{k - 1} \times N_{1}} & \cdots & A_{k - 1,k - 1} & A_{k - 1,k} \\
0_{M_{k} \times N_{1}} & \cdots & 0_{M_{k} \times N_{k - 1}} & A_{kk} \\
 & & & 
\end{pmatrix} - \begin{pmatrix}
0_{M_{1} \times N_{1}} & \cdots & B_{1,k - 1} & B_{1,k} \\
\vdots & \ddots & \vdots & \vdots \\
0_{M_{k - 1} \times N_{1}} & \cdots & 0_{M_{k - 1} \times N_{k - 1}} & B_{k - 1,k} \\
0_{M_{k} \times N_{1}} & \cdots & 0_{M_{k} \times N_{k - 1}} & B_{kk} \\
 & & & 
\end{pmatrix}z$$ The diagonal blocks are of dimension
$\left( M_{i} \times N_{i} \right)$. The first $(k - 1)$ diagonal blocks
of $A$ have full row rank (and hence $0 \leq M_{i} \leq N_{i}$) and the
first $(k - 1)$ diagonal block of $B$ are zero. Therefore $P(z)$ has
full row rank if and only if the last block
$\left( A_{kk} - B_{kk}z \right)$ has full row rank. For this last block
the following cases may occur:

1.  If $\left( A_{kk} - B_{kk} \right)$ is “tall” ($M_{k} > N_{k}$),
    then $P(z)$ has a non trivial left kernel for all
    $z \in {\mathbb{C}}$.  
    In this case $P$ (and the matrix $c$) is *not* left prime.
2.  If the pencil $\left( A_{kk} - B_{kk} \right)$ is square
    ($M_{k} = N_{k}$) and $B_{kk}$ is non singular, then the zeroes of
    the pencil $P(z)$ (and the zeroes of $c(z)$) are the eigenvalues of
    the matrix $A_{kk}B_{kk}^{- 1}$.  
    In this case $P$ (and the matrix $c$) is *not* left prime.
3.  If $B_{kk}$ is zero and the matrix $A$ has full row rank, then
    $P(z)$ has full row rank for all $z \in {\mathbb{C}}$.  
    In this case the pencil $P$ (and the matrix $c$) is left prime.
4.  If $B_{kk}$ is zero and the matrix $A$ does not have full row rank,
    then $P(z)$ has a non trivial left kernel for all
    $z \in {\mathbb{C}}$.  
    In this case $P$ (and the matrix $c$) is *not* left prime.

The function
[`is.coprime()`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md)
uses this approach to test wether a polynomial is left prime (or a pair
of polynomials is left coprime). In the default case the function just
returns a boolean variable. However, if the optional argument
`only.answer` is set to `FALSE`, then a list is returned which contains
the above described staircase form. The list contains the following
slots:

`answer` is `TRUE`, if the polynomial is prime (the pair is left
coprime).  
`A`,`B` hold the two matrices $A$ and $B$ (in staircase form).  
`m`,`n` are two integer vectors which contain the size of the blocks
$\left( M_{i} \right)$ and $\left( N_{i} \right)$ respectively.  
`zeroes` is a vector which contains the zeroes of the pencil (the
polynomial $c$). In the (co-)prime case this vector is empty and if $P$
is rank deficient for all $z \in {\mathbb{C}}$ then `zeroes = NA_real_`.

### Examples

We present three examples pertaining to the three different outcomes:

- coprime polynomial matrices $a(z)$ and $b(z)$
- non-coprime polynomial matrices with a finite number of common zeros
- non-coprime polynomial matrices with reduced rank for all
  $z \in {\mathbb{C}}$

#### Coprime Polynomials

Two “generic” polynomials are coprime.

``` r
# Generate a random (m x m), and a random (m x n) polynomial matrix, with degree p=2
m = 2
n = 3
set.seed(1803)
a = test_polm(dim = c(m,m), degree = 2, digits = 1, random = TRUE)
b = test_polm(dim = c(m,n), degree = 2, digits = 1, random = TRUE)
```

##### Hermite Form

The Hermite form of $\left( a(z),b(z) \right)$ is

``` r
HNF = hnf(cbind(a,b), from_left = FALSE)
print(HNF$h, digits = 2, format = 'character')
#> ( 2 x 5 ) matrix polynomial with degree <= 0 
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,]     1     0     0     0     0
#> [2,]     0     1     0     0     0
```

The square submatrix $h_{1}$ of the Hermite form $h$ is equal to the
identity matrix, hence $(a,b)$ are left coprime.

##### Pencil

The method using the pencil gives the same answer:

- There are no common zeros.
- The polynomial matrices are left coprime.

``` r
out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
cat("The polynomials are left coprime: ", out$answer,"\n", 
    "The zeros of the pencil are: ", out$zeroes, sep = '')
#> The polynomials are left coprime: TRUE
#> The zeros of the pencil are:
```

#### Not coprime: Finite Number of Zeros

Here we consider an example where
$c(z) = \left\lbrack a(z),b(z) \right\rbrack$ is rank deficient for some
(but not all) $z \in {\mathbb{C}}$. We generate these polynomials by
multiplying them with a common (matrix) factor. We show that (at least)
one singular value of $c(z)$ evaluated at the zeros of the common factor
is zero.

``` r
a0 = a
b0 = b

# Generate random common factor  with degree 1
r = test_polm(dim = c(m,m), degree = 1, digits = 1, random = TRUE)

# Generate polynomials with a common factor
a = r %r% a0
b = r %r% b0
c = cbind(a,b)

z_r = zeroes(r)
cat("The zeros of the common factor r(z) are: \n", z_r, "\n\n")
#> The zeros of the common factor r(z) are: 
#>  -0.2051932 5.377607

d = svd(unclass(zvalues(c, z = z_r[1]))[,,1])$d
cat("minimal singular value of c(z_0): ", d[m], sep="")
#> minimal singular value of c(z_0): 1.007532e-16
```

##### Hermite Form

The Hermite Form of $c(z) = \left( a(z),b(z) \right)$ is

``` r
HNF = hnf(cbind(a,b), from_left = FALSE)
print(HNF$h, digits = 2, format = 'character')
#> ( 2 x 5 ) matrix polynomial with degree <= 2 
#>               [,1]                [,2]  [,3]  [,4]  [,5]
#> [1,]             1                   0     0     0     0
#> [2,]  2.49 - 0.39z  -1.1 - 5.17z + z^2     0     0     0
```

We calculate the zeros of $c(z) = \left\lbrack a(z),b(z) \right\rbrack$
by calculating the zeros of the diagonal elements of the square
submatrix $h_{1}(z)$. Note that for this example the $(1,1)$ element is
equal to $1$.

The (common) zeros are therefore

``` r
zeroes(HNF$h[m,m])
#> [1] -0.2051932  5.3776070
```

Note that these zeroes are identical to the zeroes of the common factor
$r(z)$.

##### Pencil

The outcome for the pencil is the same as above.

``` r
out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
cat("The polynomials are left coprime: ", out$answer,"\n", 
    "The zeros of the pencil are: ", out$zeroes, sep = '')
#> The polynomials are left coprime: FALSE
#> The zeros of the pencil are: 5.377607-0.2051932
```

#### Not coprime: Everywhere Rank Deficient

Finally let us consider the case, where
$c(z) = \left( a(z),b(z) \right)$ is rank deficient for all $z$. We
generate a pair of polynomial matrices of this kind by multiplying the
two polynomials $a(z),b(z)$ (from left) with a singular common factor.

To verify that $c(z)$ is of reduced row rank for all
$z \in {\mathbb{C}}$, we print the singular values of $c(z)$ at a
randomly selected point.

``` r
# generate a square polynomial matrix with rank one!
r = test_polm(dim = c(m,1), degree = 1, digits = 1, random = TRUE) %r% 
    test_polm(dim = c(1,m), degree = 1, digits = 1, random = TRUE) 
print(r, format = 'c')
#> ( 2 x 2 ) matrix polynomial with degree <= 2 
#>                         [,1]                     [,2]
#> [1,]   1.36 - 1.43z - 2.1z^2  -0.72 + 2.07z - 1.35z^2
#> [2,]  1.53 + 3.47z + 1.82z^2  -0.81 - 0.36z + 1.17z^2

a = r %r% a0
b = r %r% b0

# Evaluate c(z) = (a(z),b(z)) at a random point z0 and print the corresponding singular values
z0 = complex(real = rnorm(1), imaginary = rnorm(1))
d = svd( unclass(zvalues(cbind(a,b), z = z0))[,,1], nu = 0, nv = 0 )$d

cat("The singular values of c(z) \n evaluated at z0 = ", z0, "\n are: ", d, "\n\n")
#> The singular values of c(z) 
#>  evaluated at z0 =  0.3368112-0.2077692i 
#>  are:  8.087889 6.256749e-16
```

##### Hermite Form

The Hermite Form of $c(z) = \left( a(z),b(z) \right)$ is

``` r
HNF = hnf(cbind(a,b), from_left = FALSE)
print(HNF$h, digits = 1, format = 'char')
#> ( 2 x 5 ) matrix polynomial with degree <= 1 
#>              [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,]     -0.5 + z     0     0     0     0
#> [2,]  -0.6 - 0.9z     0     0     0     0
```

The matrix $c(z)$ is rank deficient (for all $z$) since
$h\lbrack 2,2\rbrack = 0$! The procedure `hnf` returns also an estimate
of the rank of $c(z)$ (when $c(z)$ is considered as rational matrix).
Here we get

``` r
cat(HNF$rank,'\n')
#> 1
```

which again means that $c(z)$ is rank deficient for all
$z \in {\mathbb{C}}$.

##### Pencil

The result for pencils is the same.

``` r
out = is.coprime(a, b, debug = FALSE, only.answer = FALSE)
cat("The polynomials are left coprime: ", out$answer,"\n", 
    "The zeros of the pencil are: ", out$zeroes, sep = '')
#> The polynomials are left coprime: FALSE
#> The zeros of the pencil are: NA
```

## Column Reduced Matrix

### Functional Area: Polynomial Manipulation (#4)

See `R/polm_methods.R` for implementation.

**Note:** The functions `degree` (`col_end_matrix`) and `col_reduce` use
different strategies to compute the column degrees (column end matrix).

- The first two consider the elements, whereas
- `col_reduce` considers the euclidean norm of the respective columns.

Let $a(z) = a^{0} + a^{1}z + \cdots + a^{p}z^{p}$,
$a^{k} = \left( a_{ij}^{k} \right)$ and
$a_{j}^{k} = \left( a_{1j}^{k},\ldots,a_{mj}^{k} \right)\prime$. The
function `degree` sets the degree of the $j$-th column to $p_{j}$ iff
$a_{ij}^{p_{j}} > \tau$ for at least one $i$ and $a_{ij}^{k} \leq \tau$
for all $i$ and $k > p_{j}$. The column end matrix is constructed
correspondingly. The function `col_reduce` sets the degree of the $j$-th
column to $p_{j}$ iff $\parallel a_{j}^{p_{j}} \parallel > \tau$ and
$\parallel a_{j}^{k} \parallel \leq \tau$ for all $k > p_{j}$.

Therefore one may get different results!

The basic strategy to construct a column reduced matrix (for the case
that $a(z)$ is square and non singular) is as follows. First permute the
columns of $a(z)$ such that the column degrees
$p_{1} \leq p_{2} \leq \cdots \leq p_{m}$ are ordered. Now suppose that
the $k$-th column, $c_{k}$ say, of the column end matrix is linearly
dependent from the previous columns ($c_{j}$, $1 \leq j < k$),
i.e. $c_{k} = \alpha_{1}c_{1} + \cdots\alpha_{k - 1}c_{k - 1}$. Then we
substract from the $k$-th column of $a(z)$,
$z^{p_{k} - p_{j}}\alpha_{j}$ times the $j$-th column. This operation
reduces the degree of the $k$-th column by one.

This procedure is repeated until the column end matrix is regular.

See also @Wolovich1974.

The coefficients $\alpha_{j}$ are determined from an SVD of the column
end matrix.

## Realization Algorithms

### Functional Area: Realization Algorithms (#2)

See `R/as_methods.R` and related functions like
[`pseries2stsp()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2stsp.md),
[`pseries2lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2lmfd.md)
for implementation.

### Impulse Response

One of the main features of the rational functions is that the Hankel
matrix of the impulse response coefficients has finite rank:

$$H = \begin{pmatrix}
k_{1} & k_{2} & \cdots \\
k_{2} & k_{3} & \cdots \\
\vdots & \vdots & 
\end{pmatrix}$$ Note that the impulse response function is only well
defined if the rational matrix has no pole at $z = 0$.

#### LMFD

We consider a $(m,n)$-dimensional rational matrix in LMFD form
$$k(z) = a^{- 1}(z)b(z) = k_{0} + k_{1}z + k_{2}z^{2} + \cdots$$ where
$a(z) = a_{0} + a_{1}z + \cdots + a_{p}z^{p}$ and
$b(z) = b_{0} + b_{1}z + \cdots + b_{p}z^{p}$. W.l.o.g we set $p = q$.
This implies $$\begin{array}{rcl}
{a(z)\left( k(z) - k_{0} \right)} & = & {b(z) - a(z)k_{0}} \\
{\left( a_{0} + a_{1}z + \cdots + a_{p}z^{p} \right)\left( k_{1}z + k_{2}z + \cdots \right)} & = & {\left( b_{0} - a_{0}k_{0} \right) + \left( b_{1} - a_{1}k_{0} \right)z + \cdots + \left( b_{p} - a_{p}k_{0} \right)z^{p}}
\end{array}$$ and thus
$$\sum\limits_{i = 0}^{p}a_{p - i}k_{j + i} = 0{\mspace{6mu}\text{for}\mspace{6mu}}j \geq 1.$$
By the above identity it follows that the Hankel matrix $H$ has finite
rank and that the coefficients of the polynomial $a(z)$ are closely
related to the left kernel of $H$. In the following we discuss, how to
construct a unique LMFD of the rational matrix $k{()}$ via the above
identity.

In order to describe the linear dependence structure of the rows of $H$
it is convenient to use a “double” index for the rows: Let
$h(i,j) \in {\mathbb{R}}^{1 \times \infty}$ denote the $j$-th row in the
$i$-th block row of $H$, i.e. $h(i,j)$ is the
$\left( (i - 1)m + j \right)$-th row of $H$.

A selection
$\mathcal{S} = \{\left( i_{k},j_{k} \right)\,|\, k = 1,\ldots,s\}$ of
rows of the Hankel matrix is called a *nice* selection, if there are no
“holes” in the sense that $(i,j) \in \mathcal{S}$, $i > 1$ implies that
$(i - 1,j) \in \mathcal{S}$. Nice selections may be described by a
multi-index $\nu = \left( \nu_{1},\ldots,\nu_{m} \right)$, where
$\nu_{j} = \max\{ i\,|\,(i,j) \in \mathcal{S}\}$.

Suppose that $H$ has rank $s$. In general there are many different
selections of rows of $H$ which form a basis for the row space of $H$.
In the following we choose the *first* $s$ rows of $H$ which form a
basis of the row space and denote the corresponding selection with
$\mathcal{S} = \{\left( i_{k},j_{k} \right)\,|\, k = 1,\ldots,s\}$. Due
to the Hankel structure of $H$ this is a nice selection in the above
sense. The corresponding $\nu_{j}$’s are called *Kronecker indices* of
the Hankel matrix (respectively of the rational matrix $k( \cdot )$).
Note that the sum of the Kronecker indices is equal to the rank of $H$:
$\sum_{j = 1}^{m}\nu_{j} = s$.

The row $h\left( \nu_{j} + 1,j \right)$ is linearly dependent on the
previous (basis) rows and therefore
$$h\left( \nu_{j} + 1,j \right) + \sum\limits_{\substack{{(\nu_{j} + 1 - i,l)} \in \mathcal{S} \\ {(\nu_{j} + 1 - i,l)} < {(\nu_{j} + 1,j)}}}a_{i,jl}h\left( \nu_{j} + 1 - i,l \right) = 0$$
holds for suitably chosen coefficients $a_{i,jl}$. The sum runs over all
basis rows previous to the row $h\left( \nu_{j} + 1,j \right)$, i.e. the
notation $(k,l) < (i,j)$ means ($k < i$) or ($k = i$ and $l < j)$). This
equation now (uniquely) defines a polynomial matrix
$a(z) = a_{0} + a_{1}z + \cdots + a_{p}z^{p}$, with
$p = \max_{j}\nu_{j}$.  
Let $a_{i,jl}$ denote the $(jl)$-th entry of the matrix $a_{i}$. The
entries $a_{i,jl}$, for
$\left( \nu_{j} + 1 - i,l \right) \in \mathcal{S}$ and
$\left( \nu_{j} + 1 - i,l \right) < \left( \nu_{j} + 1,j \right)$ are
read off from the above equation, $a_{0,jj}$ is set to one
($a_{0,jj} = 1$) and all other entries are set to zero.

By construction $a(z)\left( k(z) - k_{0} \right)$ is a polynomial (with
degree less than or equal to $p = \max_{j}\nu_{j}$). In the last step we
set $$b(z) = a(z)\left( k(z) - k_{0} \right) + a(z)k_{0}.$$

By the above construction one gets a unique LMFD representation of the
rational matrix $k(z)$ with the following properties:

- $a_{0}$ is a lower triangular matrix with ones on the diagonal.
- $b_{0} = a_{0}k_{0}$. In particular for a square rational matrix with
  $k_{0} = I_{m}$ we have $b_{0} = a_{0}$.
- the row degrees of $\left( a(z),b(z) \right)$ are equal to
  $\nu_{1},\ldots,\nu_{m}$,  
- the elements $a_{ij}(z)$ are divisible by $z^{\nu_{ij}}$, where
  $\nu_{ij} = \max\left( \nu_{i} + 1 - \nu_{j},1 \right)$ for $j > i$
  and $\nu_{ij} = \max\left( \nu_{i} + 1 - \nu_{j},0 \right)$ for
  $j < i$.
- the pair $\left( a(z),b(z) \right)$ is left coprime and row reduced.

A pair $\left( a(z),b(z) \right)$ which satisfied the above conditions
is said to be in *echelon canonical form*.

#### Ho-Kalman

A quite analogous strategy may be used to construct a (unique) state
space realization of a rational matrix. Suppose that the Hankel matrix
$H$ has rank $s$ and that $S \in {\mathbb{R}}^{s \times \infty}$ is such
that $SH$ is a basis for the row space of $H$. Then a state space
representation of the rational matrix $k( \cdot )$ is obtained by a
solving the following equations $$\begin{array}{rcl}
{SH_{- 1,.}} & = & {ASH} \\
H_{1,.} & = & {CSH} \\
{SH_{.,1}} & = & B \\
k_{0} & = & {D,}
\end{array}$$ where $$H_{- 1,.} = \begin{pmatrix}
k_{2} & k_{3} & \cdots \\
k_{3} & k_{4} & \cdots \\
\vdots & \vdots & 
\end{pmatrix}\,\quad H_{1,.} = \begin{pmatrix}
k_{1} & k_{2} & \cdots
\end{pmatrix}\,{\mspace{6mu}\text{and}\mspace{6mu}}H_{,1} = \begin{pmatrix}
k_{1} \\
k_{2} \\
\vdots
\end{pmatrix}$$

The obtained result clearly depends on the matrix $S$, i.e. on the
choice of the basis of the row space of $H$.

Two choices are implemented:

- echelon form: here $S$ is the selection matrix which selects the first
  $s$ linearly independent rows of $H$. The statespace realization
  obtained then is in *echelon canonical form*.
- balanced form: here $S$ is obtained via an SVD decomposition of a
  finite dimensional sub matrix of $H$. For more details, see …

#### Notes:

The above constructions for a LMFD or statespace representations also
work for a finite dimensional sub matrix $$H_{fp} = \begin{pmatrix}
k_{1} & k_{2} & \cdots & k_{p} \\
k_{2} & k_{3} & \cdots & k_{p + 1} \\
\vdots & \vdots & & \vdots \\
k_{f} & k_{f + 1} & \cdots & k_{f + p - 1}
\end{pmatrix}$$ provided that $(f,p)$ is “large enough”.

For a more detailed discussion on Kronecker indices and echelon
canonical forms, see @Hannan.Deistler12.

## Operations with State-Space Realizations

### Functional Area: State-Space Tools (#6)

See `R/arithmetic_methods.R` and `R/stsp_methods.R` for implementation.

### Addition

$$\left( C_{1}\left( z^{- 1} - A_{1} \right)^{- 1}B_{1} + D_{1} \right) + \left( C_{2}\left( z^{- 1} - A_{2} \right)^{- 1}B_{2} + D_{2} \right)$$

\$\$ \left\[\begin{array}{@{}cc\|c@{}} A_1 & 0 & B_1 \\ 0 & A_2 & B_2 \\
\hline C_1 & C_2 & D_1 + D_2 \end{array}\right\] \$\$

See
[`Ops.ratm()`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
in particular `a + b`.

### Multiplication

$$\left( C_{1}\left( z^{- 1} - A_{1} \right)^{- 1}B_{1} + D_{1} \right) \cdot \left( C_{2}\left( z^{- 1} - A_{2} \right)^{- 1}B_{2} + D_{2} \right)$$

\$\$ \left\[\begin{array}{@{}cc\|c@{}} A_1 & B_1 C_2 & B_1 D_2 \\ 0 &
A_2 & B_2 \\ \hline C_1 & D_1 C_2 & D_1 D_2 \end{array}\right\] \$\$

See `a %r% b`.

### Inverse

$$\left( C\left( z^{- 1} - A \right)^{- 1}B + D \right)^{- 1}$$

\$\$ \left\[\begin{array}{@{}c\|c@{}} A - B D^{-1} C & BD^{-1} \\ \hline
-D^{-1}C & D^{-1} \end{array}\right\] \$\$

See
[`Ops.ratm()`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
in particular `x^{-1}`.

### Bind by columns or rows

$$\begin{pmatrix}
{C_{1}\left( z^{- 1} - A_{1} \right)^{- 1}B_{1} + D_{1})} & {C_{2}\left( z^{- 1} - A_{2} \right)^{- 1}B_{2} + D_{2})}
\end{pmatrix}$$

\$\$ \left\[\begin{array}{@{}cc\|cc@{}} A_1 & 0 & B_1 & 0 \\ 0 & A_2 & 0
& B_2 \\ \hline C_1 & C_2 & D_1 & D_2 \end{array}\right\] \$\$ See
[`cbind.ratm()`](https://bfunovits.github.io/rationalmatrices/reference/bind.md).

$$\begin{pmatrix}
{C_{1}\left( z^{- 1} - A_{1} \right)^{- 1}B_{1} + D_{1}} \\
{C_{2}\left( z^{- 1} - A_{2} \right)^{- 1}B_{2} + D_{2}}
\end{pmatrix}$$

\$\$ \left\[\begin{array}{@{}cc\|c@{}} A_1 & 0 & B_1 \\ 0 & A_2 & B_2 \\
\hline C_1 & 0 & D_1 \\ 0 & C_2 & D_2 \end{array}\right\] \$\$ See
[`rbind.ratm()`](https://bfunovits.github.io/rationalmatrices/reference/bind.md).

### Elementwise Multiplication

First consider the product of two $(m \times 1)$ matrices, and write the
elementwise product as the product of a suitable diagonal matrix with
the second factor

\$\$ \left\[\begin{array}{@{}cccc\|c@{}} A_1 & \cdots & 0 & B_1 C\_{2,1}
& B_1 D\_{2,1} \\ \vdots & \ddots & \vdots & \vdots & \vdots \\ 0 &
\cdots & A_1 & B_1 C\_{2,m} & B_1 D\_{2,m} \\ 0 & \cdots & 0 & A_2 & B_2
\\ \hline C\_{1,1} & \cdots & 0 & D\_{1,1} C\_{2,1} & D\_{1,1} C\_{2,1}
\\ \vdots & \ddots & \vdots & \vdots & \vdots \\ 0 & \cdots & C\_{1,m} &
D\_{1,m} C\_{2,m} & D\_{1,m} D\_{2,m} \end{array}\right\] \$\$ where
$C_{i,k}$ denotes the $k$-th row of $C_{i}$ and $D_{i,k}$ is the
$k$-entry of the $(m \times 1)$-dimensional matrix $D_{i}$.

It can be shown (???) that the controllability matrix of the above
statespace model has rank at most $\left( s_{1} + s_{2} \right)$.
Therefore we may construct an equivalent model with state dimension
$\left( s_{1} + s_{2} \right)$.

Now we simply do this construction for all columns of the two factors
$a(z)$ and $b(z)$ and then “column-bind” the result. The statespace
realisation then has state dimension $n\left( s_{1} + s_{2} \right)$.

For “wide” matrices ($m < n$) we consider the transpose of the two
factors. This gives statespace realisations with
$$s = \min(m,n)\left( s_{1} + s_{2} \right).$$

- proof ???
- kann man die $\left( s_{1} + s_{2} \right)$ Zustandsraum Darstellung
  direkt “hinschreiben”?

See
[`Ops.ratm()`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
in particular `a * b`.

### Transpose

$$\left( C\left( z^{- 1}I - A \right)^{- 1}B + D \right)\prime = B\prime\left( z^{- 1}I - A\prime \right)^{- 1}C + D\prime$$

\$\$ \left\[\begin{array}{@{}c\|c@{}} A' & C' \\ \hline B' & D'
\end{array}\right\] \$\$

See [`t()`](https://rdrr.io/r/base/t.html).

### Hermitean Transpose

$$\begin{array}{rcl}
\left\lbrack C\left( z^{- 1}I - A \right)^{- 1}B + D \right\rbrack^{*} & = & {B\prime(zI - A\prime)^{- 1}C\prime + D\prime} \\
 & = & {B\prime\left\lbrack ( - A\prime z)\left( z^{- 1}I - A^{- T} \right) \right\rbrack^{- 1}C\prime + D\prime} \\
 & = & {B\prime\left\lbrack \left( - A^{- T}z^{- 1} \right)\left( zI + z^{2}A^{- T} + z^{3}A^{- 2T} + \cdots \right) \right\rbrack C\prime + D\prime} \\
 & = & {\left( - A^{- 1}B \right)\prime\left\lbrack zI + z^{2}A^{- T} + Z^{3}A^{- 2T} + \cdots) \right\rbrack\left( CA^{- 1} \right)\prime + D\prime - B\prime A^{- T}C\prime} \\
 & = & {\left( - A^{- 1}B \right)\prime\left( z^{- 1}I - A^{- T} \right)^{- 1}\left( CA^{- 1} \right)\prime + \left( D\prime - B\prime A^{- T}C\prime \right)} \\
 & & 
\end{array}$$

\$\$ \left\[\begin{array}{@{}c\|c@{}} A^{-T} & (CA^{-1})' \\ \hline
(-A^{-1}B)' & D' - B' A^{-T} C' \end{array}\right\] \$\$

See
[`Ht()`](https://bfunovits.github.io/rationalmatrices/reference/Ht.md).

### Derivative

$$\begin{array}{rcl}
{\frac{d}{dz}\left( C\left( z^{- 1} - A \right)^{- 1}B + D \right)} & = & {\frac{d}{dz}\left( D + CBz + CABz^{2} + CA^{2}Bz^{3} + CA^{2}Bz^{4} + \cdots \right)} \\
 & = & {CB + 2CABz + 3CA^{2}Bz^{2} + 4CA^{2}Bz^{3} + \cdots} \\
 & = & {\bar{C}\bar{A}\bar{B} + \bar{C}{\bar{A}}^{2}\bar{B}z + \bar{C}{\bar{A}}^{3}\bar{B}z^{2} + \bar{C}{\bar{A}}^{4}\bar{B}z^{3} + \cdots} \\
 & & 
\end{array}$$ with $$\bar{C} = \begin{pmatrix}
C & 0
\end{pmatrix},\ \bar{A} = \begin{pmatrix}
A & I \\
0 & A
\end{pmatrix}{\mspace{6mu}\text{and}\mspace{6mu}}\bar{B} = \begin{pmatrix}
0 \\
B
\end{pmatrix}.$$ The corresponding statespace realization is

\$\$ \left\[\begin{array}{@{}cc\|c@{}} A & I & B \\ 0 & A & AB \\ \hline
CA & C & CB \end{array}\right\] \$\$ See
[`derivative.stsp()`](https://bfunovits.github.io/rationalmatrices/reference/derivative.md).

### Construct statespace representation of rational matrix in LMFD form

Suppose $c(z) = a^{- 1}(z)b(z)$ with
$a(z) = a_{0} + a_{1}z + a_{2}z^{2} + \cdots + a_{p}z^{p}$,
$b(z) = b_{0} + b_{1}z + \cdots + b_{q}z^{q}$ is given. The
construction, detailed below, assumes that $a_{0}$ is non singular and
hence we rewrite the polynomial $a(z)$ as  
$a(z) = I_{m} - a_{1}z - a_{2}z^{2} - \cdots - a_{p}z^{p}$. Furthermore
we assume, w.l.o.g that $p = q$. The powerseries expansion of
$c( \cdot )$ is $c(z) = k_{0} + k_{1}z + \cdots + k_{p}z^{p} + \cdots$.

Then a statespace realization of $c( \cdot )$ is \$\$
\left\[\begin{array}{@{}cccc\|c@{}} a_1 & \cdots & a\_{p-1} & a_p & k_p
\\ I_m & \cdots & 0 & 0 & k\_{p-1} \\ \vdots& \ddots & \vdots & \vdots &
\vdots \\ 0 & \cdots & I_m & 0 & k_1 \\ \hline 0 & \cdots & 0 & I_m &
k_0 \end{array}\right\] \$\$ This scheme is implemented in
`as.stsp.lmfd` (and `as.stsp.polm`).

### Construct statespace representation of rational matrix in RMFD form

Likewise, it is possible to write an RMFD \$k(z)=d(z) c^{-1}(z) \$,
where $c(z) = c_{0} + c_{1}z + c_{2}z^{2} + \cdots + c_{p}z^{p}$, and
$d(z) = d_{0} + d_{1}z + \cdots + d_{q}z^{q}$ is given. Again, we
simplify the exposition by setting the non-singular $c_{0}$ equal to the
$m$-dimensional identity matrix $I_{m}$. In the case $p \geq q$, the
state space realization of $k(z)$ is

\$\$ \left\[\begin{array}{@{}cccc\|c@{}} c_1 & \cdots & c\_{p-1} & c_p &
I_m \\ I_m & \cdots & 0 & 0 & 0 \\ \vdots& \ddots & \vdots & \vdots &
\vdots \\ 0 & \cdots & I_m & 0 & 0 \\ \hline d_1 & \cdots & d\_{q-1} &
d_q & d_0 \end{array}\right\] \$\$ If $p > q$, some zeros appear in the
$C$ matrix. If $p < q$, some zeros are padded into the $A$ matrix. This
scheme is implemented in `as.stsp.rmfd`.

## Schur, QZ Decomposition, Lyapunov and Riccati Equations

### Functional Area: Numerical Methods (#7)

See `src/lyapunov.cpp`, `inst/include/rationalmatrices_lyapunov.h`, and
`R/lyapunov.R` for implementation.

### Lyapunov Equation

The generalized (non-symmetric) Lyapunov equation is

$$X = AXB^{*} + Q$$

With a Schur factorization of $A$ and $B$ we obtain a form with
triangular matrices $A$, $B$

$$\begin{array}{rcl}
\begin{pmatrix}
X_{11} & X_{12} \\
X_{21} & X_{22}
\end{pmatrix} & = & {\begin{pmatrix}
A_{11} & A_{12} \\
0 & A_{22}
\end{pmatrix}\begin{pmatrix}
X_{11} & X_{12} \\
X_{21} & X_{22}
\end{pmatrix}\begin{pmatrix}
B_{11}^{*} & 0 \\
B_{12}^{*} & B_{22}^{*}
\end{pmatrix} + \begin{pmatrix}
Q_{11} & Q_{12} \\
Q_{21} & Q_{22}
\end{pmatrix}} \\
 & = & {\begin{pmatrix}
{A_{11}X_{11}B_{11}^{*} + A_{11}X_{12}B_{12}^{*} + A_{12}X_{21}B_{11}^{*} + A_{12}X_{22}B_{12}^{*}} & {\left( A_{11}X_{12} + A_{12}X_{22} \right)B_{22}^{*}} \\
{A_{22}\left( X_{21}B_{11}^{*} + X_{22}B_{12}^{*} \right)} & {A_{22}X_{22}B_{22}^{*}}
\end{pmatrix} + \begin{pmatrix}
Q_{11} & Q_{12} \\
Q_{21} & Q_{22}
\end{pmatrix}}
\end{array}$$

- 22-block: solve for $X_{22}$
- 12-block: solve for $X_{12}$
- 21-block: solve for $X_{21}$
- continue with 11-block

we could also consider the non-square case, where $A$ is $m \times m$,
$B$ is $n \times n$ and $Q$ is $m \times n$.

### Balancing and Balanced Truncation

Consider a pair of Grammians $P,Q$. We first compute the sqare roots
$P = MM\prime$ and $Q = NN\prime$ and then consider the SVD[¹](#fn1) of
$M\prime N$:

$$M\prime N = U\Sigma V\prime = \left\lbrack U_{1},U_{2} \right\rbrack\begin{pmatrix}
\Sigma_{11} & 0 \\
0 & \Sigma_{22}
\end{pmatrix}\left\lbrack V_{1},V_{2} \right\rbrack\prime$$ where the
diagonal blocks $\Sigma_{11}$, $\Sigma_{22}$ are of size
$\left( s_{1},s_{1} \right)$ and $\left( s_{2},s_{2} \right)$
respectively. We assume that $\Sigma_{11}$ is positive definite and set
$$\begin{array}{rcl}
T_{1} & = & {\Sigma_{11}^{- 1/2}V_{1}\prime N\prime} \\
S_{1} & = & {MU_{1}\Sigma_{11}^{- 1/2}}
\end{array}$$ and note that $$\begin{array}{rcl}
{T_{1}S_{1}} & = & {\Sigma_{11}^{- 1/2}V_{1}\prime N\prime MU_{1}\Sigma_{11}^{- 1/2}} \\
 & = & {\Sigma_{11}^{- 1/2}V_{1}\prime V\Sigma U\prime MU_{1}\Sigma_{11}^{- 1/2} = I_{s_{1}}}
\end{array}$$ We extend these matrices to square $(s \times s)$ matrices
$$T = \begin{pmatrix}
T_{1} \\
T_{2}
\end{pmatrix}{\mspace{6mu}\text{and}\mspace{6mu}}S = \begin{pmatrix}
S_{1} & S_{2}
\end{pmatrix}$$

such that $TS = ST = I_{s}$, i.e. $S = T^{- 1}$. To this end we consider
the SVD’s

$$\begin{array}{rclcl}
{S_{1}T_{1}} & = & {{\bar{U}}_{1}{\bar{\Sigma}}_{11}{\bar{V}}_{1}\prime} & & {{\bar{U}}_{1},{\bar{V}}_{1} \in {\mathbb{R}}^{s \times s_{1}},\,{\bar{\Sigma}}_{11} \in {\mathbb{R}}^{s_{1} \times s_{1}}} \\
{{\bar{U}}_{2}\prime{\bar{V}}_{2}} & = & {\widehat{U}\widehat{\Sigma}\widehat{V}\prime} & & {\widehat{U},\widehat{V},\widehat{\Sigma} \in {\mathbb{R}}^{s_{2} \times s_{2}}}
\end{array}$$ and set
$T_{2} = {\widehat{\Sigma}}^{- 1/2}\widehat{U}\prime{\bar{U}}_{2}\prime$
and $S_{2} = {\bar{V}}_{2}\widehat{V}{\widehat{\Sigma}}^{- 1/2}$.

We now consider the *transformed* statespace realization \$\$
\left\[\begin{array}{@{}c\|c@{}} A & B \\ \hline C & D
\end{array}\right\] \\ \longrightarrow \\
\left\[\begin{array}{@{}cc\|c@{}} T_1 A S_1 & T_1 A S_2 & T_1 B \\ T_2 A
S_1 & T_2 A S_2 & T_2 B \\ \hline C S_1 & C S_2 & D \end{array}\right\]
\$\$ and the Grammians are transformed as
$\left. P\rightarrow TPT\prime \right.$ and
$\left. Q\rightarrow S\prime QS \right.$.

$$TM = \begin{pmatrix}
{\Sigma_{11}^{- 1/2}V_{1}\prime N\prime M} \\
{T_{2}M}
\end{pmatrix} = \begin{pmatrix}
{\Sigma_{11}^{1/2}U_{1}\prime} \\
{T_{2}M}
\end{pmatrix}$$

$$S\prime N = \begin{pmatrix}
{\Sigma_{11}^{- 1/2}U_{1}\prime M\prime N} \\
{S_{2}\prime N}
\end{pmatrix} = \begin{pmatrix}
{\Sigma_{11}^{1/2}V_{1}\prime} \\
{S_{2}\prime N}
\end{pmatrix}$$$$TPT\prime = (TM)(TM)\prime = \begin{pmatrix}
{\Sigma_{11}^{1/2}U_{1}\prime U_{1}\Sigma_{11}^{1/2}} & {\Sigma_{11}^{1/2}U_{1}\prime M\prime T_{2}\prime} \\
{T_{2}MU_{1}\Sigma_{11}^{1/2}} & {T_{2}MM\prime T_{2}\prime}
\end{pmatrix} = \begin{pmatrix}
\Sigma_{11} & 0 \\
0 & {T_{2}PT_{2}\prime}
\end{pmatrix}$$

$$S\prime QS = (S\prime N)(S\prime N)\prime = \begin{pmatrix}
{\Sigma_{11}^{- 1/2}V_{1}\prime V_{1}\Sigma_{11}^{1/2}} & {\Sigma_{11}^{1/2}V_{1}\prime N\prime S_{2}} \\
{S_{2}\prime NV_{1}\Sigma_{11}^{1/2}} & {S_{2}\prime NN\prime S_{2}}
\end{pmatrix} = \begin{pmatrix}
\Sigma_{11} & 0 \\
0 & {S_{2}\prime QS_{2}}
\end{pmatrix}$$ Here we have used that
$T_{2}MU_{1}\Sigma_{11}^{1/2} = T_{2}S_{1}\Sigma_{11} = 0$ and
$\Sigma_{11}^{1/2}V_{1}\prime N\prime S_{2} = \Sigma_{11}T_{1}S_{2} = 0$.
Due to the block diagonal structure of $TPT\prime$ and $S\prime QS$ we
also immediately get

$$TPQS = TPT\prime S\prime QS = \begin{pmatrix}
\Sigma_{11}^{2} & 0 \\
0 & {T_{2}PQS_{2}}
\end{pmatrix}$$

and

$$T_{2}PQS_{2} = T_{2}MM\prime NN\prime S_{2} = T_{2}MU\Sigma V\prime N\prime S_{2} = \underset{= 0}{\underbrace{T_{2}MU_{1}}}\Sigma_{11}\underset{= 0}{\underbrace{V_{1}\prime N\prime S_{2}}} + T_{2}MU_{2}\Sigma_{22}V_{2}\prime N\prime S_{2} = T_{2}MU_{2}\Sigma_{22}V_{2}\prime N\prime S_{2}.$$

The following scenarios are of particular interest:

- $s_{2} = 0$: The model is *minimal* and the above procedure renders
  the statespace realization into balanced form.
- $s_{2} > 0$ and $\Sigma_{22} = 0$: The above procedure renders the
  model in a form where the controllable and observable states are
  clearly seperated from the non observable or non controllable states.
  In particular the truncated model \$\$
  \left\[\begin{array}{@{}c\|c@{}} T_1 A S_1 & T_1 B \\ T_2 A S_1 & T_2
  B \\ \hline C S_1 & D \end{array}\right\] \$\$ is equivalent to the
  original model and it is minimal and in balanced form
- If $s_{2} > 0$ but $\Sigma_{22} \neq 0$ then the truncated model is
  just an approximation of the original model. The quality depends on
  the size of the neglected singular values. Note also that the
  truncated model is not balanced.

Note that the determination of the number of non zeroe singular values
is quite tricky. (difficult due to numerics). Currently the number of
non zero singular values is chosen as the number of singular values
which satisfy $\sigma_{k} > \tau\sigma_{1}$ where $\tau$ is a user
defined tolerance parameter.

### Reflect Poles and Zeroes

#### Functional Area: Pole/Zero Reflection (#8)

See `R/reflect_poles_zeroes.R` for implementation.

#### All-Pass Matrices

A rational matrix $k(z)$ is all-pass if
$k(z)k\prime\left( \frac{1}{z} \right) = k\prime\left( \frac{1}{z} \right)k(z) = I_{m}$.
If $(A,B,C,D)$ is a statespace realization of $k(z)$ then the product
$k(z)k\prime\left( \frac{1}{z} \right)$ has a realization given by \$\$
\left\[\begin{array}{@{}cc\|c@{}} A & -BB'A^{-T} & BD' - BB'A^{-T}C' \\
0 & A^{-T} & A^{-T}C' \\ \hline C & -DB'A^{-T} & DD' - DB' A^{-T} C'
\end{array}\right\] \$\$

A state transformation gives

\$\$ \left\[ \begin{array}{@{}cc\|c@{}} I & X & 0 \\ 0 & I & 0 \\ \hline
0 & 0 & I \end{array} \right\] \left\[\begin{array}{@{}cc\|c@{}} A &
-BB'A^{-T} & BD' - BB'A^{-T}C' \\ 0 & A^{-T} & A^{-T}C' \\ \hline C &
-DB'A^{-T} & DD' - DB' A^{-T} C' \end{array}\right\] \left\[
\begin{array}{@{}cc\|c@{}} I & -X & 0 \\ 0 & I & 0 \\ \hline 0 & 0 & I
\end{array} \right\] = \left\[\begin{array}{@{}cc\|c@{}} A &
XA^{-T}-AX-BB'A^{-T} & BD' - BB'A^{-T}C' + X A^{-T}C'\\ 0 & A^{-T} &
A^{-T}C' \\ \hline C & -DB'A^{-T}-CX & DD' - DB' A^{-T} C'
\end{array}\right\] \$\$

The (1,1) block is not controllable and the (2,2) block is not
observable if

$$\begin{array}{rclcrcl}
{XA^{- T} - AX - BB\prime A^{- T}} & = & 0 & \Leftrightarrow & {X - AXA\prime - BB\prime} & = & 0 \\
{- DB\prime A^{- T} - CX} & = & 0 & \Leftrightarrow & {\left\lbrack X\prime,A^{- 1}B \right\rbrack\lbrack C,D\rbrack\prime} & = & 0 \\
{BD\prime - BB\prime A^{- T}C\prime + XA^{- T}C\prime} & = & 0 & \Leftrightarrow & {\left\lbrack XA^{- T} - BB\prime A^{- T},B \right\rbrack\lbrack C,D\rbrack\prime} & = & 0 \\
 & & & \Leftrightarrow & {\lbrack AX\prime,B\rbrack\lbrack C,D\rbrack\prime} & = & 0 \\
 & & & & & & 
\end{array}$$ Furthermore
$$DD\prime - DB\prime A^{- T}C\prime = D\left\lbrack B\prime A^{- T}, - I \right\rbrack\lbrack C,D\rbrack\prime = I$$
must hold.

It follows that we end up with the system \$\$
\left\[\begin{array}{@{}cc\|c@{}} A & 0 & 0\\ 0 & A^{-T} & A^{-T}C' \\
\hline C & 0 & I \end{array}\right\] \$\$ whose transfer function is
equal to $$\begin{pmatrix}
C & 0
\end{pmatrix}\left( Iz^{- 1} - \begin{pmatrix}
A & 0 \\
0 & A^{- T}
\end{pmatrix} \right)^{- 1}\begin{pmatrix}
0 \\
{A^{- T}C\prime}
\end{pmatrix} + I_{m},$$ i.e. the identity matrix as required.

In the following we need to construct an allpass matrix with given
$(A,B)$ (or given $(A,C)$).

1.  From the above Lyapunov equation $X = AXA\prime + BB\prime$
    determine $X$.
2.  Determine the row space of $\lbrack C,D\rbrack$ from the left kernel
    of $\lbrack AX\prime,B\rbrack\prime$.
3.  Determine the scaling matrix from the requirement
    $DD\prime - DB\prime A^{- T}C\prime = I$.

*Note:* Wann gibt es Probleme? ….

#### zero cancellations

Let $k(z)$ be given. We want to find an allpass transfer function
$\widehat{k}(z)$ such some of the poles of $k(z)$ are cancelled in
$k(z)\widehat{k}(z)$.

Consider the state transformation of the realization of $k\widehat{k}$:

\$\$ \left\[ \begin{array}{@{}ccc\|c@{}} I & 0 & 0 & 0 \\ 0 & I & 0 & 0
\\ I & 0 & I & 0 \\ \hline 0 & 0 & 0 & I \end{array} \right\] \left\[
\begin{array}{@{}ccc\|c@{}} A\_{11} & A\_{12} & B\_{1}\hat{C} &
B\_{1}\hat{D} \\ A\_{21} & A\_{22} & B\_{2}\hat{C} & B\_{2}\hat{D} \\ 0
& 0 & \hat{A} & \hat{B} \\ \hline C\_{1} & C\_{2} & D\hat{C} & D\hat{D}
\end{array} \right\] \left\[ \begin{array}{@{}ccc\|c@{}} I & 0 & 0 & 0
\\ 0 & I & 0 & 0 \\ -I & 0 & I & 0 \\ \hline 0 & 0 & 0 & I \end{array}
\right\] = \left\[ \begin{array}{@{}ccc\|c@{}} A\_{11}-B\_{1}\hat{C} &
A\_{12} & B\_{1}\hat{C} & B\_{1}\hat{D}\\ A\_{21}-B\_{2}\hat{C} &
A\_{22} & B\_{2}\hat{C} & B_2 \hat{D} \\ A\_{11}-B\_{1}\hat{C} -\hat{A}
& A\_{12} & \hat{A} + B\_{1}\hat{C} & B\_{1}\hat{D} + \hat{B} \\ \hline
C\_{1}-D\hat{C} & C\_{2} & D\hat{C} & D\hat{D} \end{array} \right\] \$\$
The $(1,1)$ block is not observable if $$\begin{array}{rclcrcl}
{A_{21} - B_{2}\widehat{C}} & = & 0 & \Leftrightarrow & {A_{21} - B_{2}D^{- 1}C_{1}} & = & 0 \\
{A_{11} - B_{1}\widehat{C} - \widehat{A}} & = & 0 & \Leftrightarrow & \widehat{A} & = & {A_{11} - B_{1}D^{- 1}C_{1}} \\
{C_{1} - D\widehat{C}} & = & 0 & \Leftrightarrow & \widehat{C} & = & {D^{- 1}C_{1}} \\
 & & & & & & 
\end{array}$$

#### pole cancellations

Let $k(z)$ be given. We want to find an allpass transfer function
$\widehat{k}(z)$ such some of the poles of $k(z)$ are cancelled in
$k(z)\widehat{k}(z)$.

Consider the state transformation of the realization of $k\widehat{k}$:

\$\$ \left\[ \begin{array}{@{}ccc\|c@{}} I & 0 & -I & 0 \\ 0 & I & 0 & 0
\\ 0 & 0 & I & 0 \\ \hline 0 & 0 & 0 & I \end{array} \right\] \left\[
\begin{array}{@{}ccc\|c@{}} A\_{11} & A\_{12} & B\_{1}\hat{C} &
B\_{1}\hat{D} \\ A\_{21} & A\_{22} & B\_{2}\hat{C} & B\_{2}\hat{D} \\ 0
& 0 & \hat{A} & \hat{B} \\ \hline C\_{1} & C\_{2} & D\hat{C} & D\hat{D}
\end{array} \right\] \left\[ \begin{array}{@{}ccc\|c@{}} I & 0 & I & 0
\\ 0 & I & 0 & 0 \\ 0 & 0 & I & 0 \\ \hline 0 & 0 & 0 & I \end{array}
\right\] = \left\[ \begin{array}{@{}ccc\|c@{}} A\_{11} & A\_{12} &
A\_{11} + B\_{1}\hat{C} - \hat{A} & B\_{1}\hat{D} - \hat{B} \\ A\_{21} &
A\_{22} & A\_{21} + B\_{2}\hat{C} & B_2 \hat{D} \\ 0 & 0 & \hat{A} &
\hat{B} \\ \hline C\_{1} & C\_{2} & C_1 + D\hat{C} & D\hat{D}
\end{array} \right\] \$\$

The (1,1) block is not controllable if

$$\begin{array}{rclcrcl}
A_{12} & = & 0 & \Leftrightarrow & A_{12} & = & 0 \\
\widehat{A} & = & {A_{11} + B_{1}\widehat{C}} & \Leftrightarrow & \left( \widehat{A} - \widehat{B}{\widehat{D}}^{- 1}\widehat{C} \right) & = & A_{11} \\
\widehat{B} & = & {B_{1}\widehat{D}} & \Leftrightarrow & \left( \widehat{B}{\widehat{D}}^{- 1} \right) & = & B_{1}
\end{array}$$

Hence the “$A,B$” matrices of the inverse of ${\widehat{k}}^{- 1}$ are
determined.

## QR Decomposition, Rank and Null Space

### Functional Area: Utilities & Helpers (#10)

Supporting theory for rank computation and null space analysis in
polynomial operations.

Consider the following example

``` r
tol = 1e-7
tol2 = tol/10

x = diag(c(1,tol2, tol2^2, tol2^3))
qr_x = qr(x, tol = tol)

x
#>      [,1]  [,2]  [,3]  [,4]
#> [1,]    1 0e+00 0e+00 0e+00
#> [2,]    0 1e-08 0e+00 0e+00
#> [3,]    0 0e+00 1e-16 0e+00
#> [4,]    0 0e+00 0e+00 1e-24
qr_x$rank
#> [1] 4
qr_x$pivot
#> [1] 1 2 3 4
qr.R(qr_x)
#>      [,1]   [,2]   [,3]  [,4]
#> [1,]   -1  0e+00  0e+00 0e+00
#> [2,]    0 -1e-08  0e+00 0e+00
#> [3,]    0  0e+00 -1e-16 0e+00
#> [4,]    0  0e+00  0e+00 1e-24
```

Is this the desired result for the rank (and the null space)? The main
reason for this somewhat strange results is that
[`qr()`](https://rdrr.io/r/base/qr.html) considers a column as zero, if
and only if the projection onto the space spanned by the previous
columns reduces the norm by a factor which is larger than `tol`.

The following example works as desired

``` r
x[1,] = 1
x[nrow(x), ncol(x)] = 1
qr_x = qr(x, tol = tol)

x
#>      [,1]  [,2]  [,3] [,4]
#> [1,]    1 1e+00 1e+00    1
#> [2,]    0 1e-08 0e+00    0
#> [3,]    0 0e+00 1e-16    0
#> [4,]    0 0e+00 0e+00    1
qr_x$rank
#> [1] 2
qr_x$pivot
#> [1] 1 4 2 3
qr.R(qr_x)
#>      [,1] [,2]   [,3]   [,4]
#> [1,]   -1   -1 -1e+00 -1e+00
#> [2,]    0   -1  0e+00  0e+00
#> [3,]    0    0 -1e-08  0e+00
#> [4,]    0    0  0e+00  1e-16
```

## References

------------------------------------------------------------------------

1.  The squared singular values of $M\prime N$ are the eigenvalues of
    $M\prime NN\prime M$. Therefore the non zero, squared singular
    values are equal to the non zero eigenvalues of
    $MM\prime NN\prime = PQ$.
