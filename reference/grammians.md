# Grammians

The procedure computes "grammians" of a statespace realization, which
may e.g. be used for balancing the statespace realization.

## Usage

``` r
grammians(
  obj,
  which = c("lyapunov", "minimum phase", "ctr", "obs", "obs_inv", "ctr_inv")
)
```

## Arguments

- obj:

  ([`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object) rational matrix in statespace form.

- which:

  (character string) specifies the type of Grammian(s) to be
  computed.See below for more details.

## Value

Either the selected Grammian (if `which` is one of
`'ctr', 'obs', 'ctr_inv', 'obs_inv'`) or a list with two components `P`
and `Q` (for the case `which = 'lyapunov'` or
`which = 'miniumum phase'`).

## Details

The *controllability Grammian* \\P\\ of a (stable) statespace
realization \$\$K(z) = C(Iz^{-1} - A)^{-1}B + D\$\$ is the solution of
the Lyapunov equation \\P = APA' + BB'\\. The *observability Grammian*
is the solution of the Lyapunov equation \\Q = A'QA + C'C\\. If the
statespace realization is *stable* (the moduli of the eigenvalues of `A`
are less than one) then \\P,Q\\ are positive semidefinite and \\P\\ is
non singular if and only if the statespace realization is *controllable*
and \\Q\\ is non singular if and only if the statespace realization is
*observable*. Hence the grammians may also be used to check whether the
statespace realization is minimal (controllable *and* observable).

If the rational matrix is (strictly) minimum phase (i.e. \\K(z)\\ is a
square, invertible matrix and the eigenvalues of the matrix \\(A -
BD^{-1}C)\\ have moduli less than one) then we may also compute the
controllability and the observability Grammian of the statespace
realization \$\$K^{-1}(z) = -D^{-1}C (Iz^{-1} - (A -
BD^{-1}C))^{-1}BD^{-1} + D^{-1}.\$\$ of the inverse matrix
\\K^{-1}(z)\\. These grammians have a similar interpretation.

The above described grammians may be selected by setting the parameter
`which` to `'ctr'`, `'obs'`, `'ctr_in'` or `'obs_inv'` respectively.

For *balancing* a statespace realization one needs a suitable pair of
grammians. Two popular choices have been implemented: For
`which = 'lyapunov'` the procedure returns the controllability and the
observability Grammian and for `which = 'minimum phase'` the
controllability matrix of the system and the observability Grammian of
the inverse system are returned.

The procedure throws an error if the state space realization is not
stable, respectively not minimum phase.

## See also

[`ctr_matrix`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md),
[`obs_matrix`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md),
[`lyapunov`](https://bfunovits.github.io/rationalmatrices/reference/lyapunov.md)
and
[`balance`](https://bfunovits.github.io/rationalmatrices/reference/balance.md).

## Examples

``` r
# create a random, (3 by 2) rational matrix, 
# with a stable and minimum phase statespüace realization
obj = test_stsp(dim = c(3,2), s = 5, bpoles = 1, bzeroes = 1)
gr = grammians(obj, which = 'lyapunov')
gr
#> $P
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,]  1.97917760  0.26060260  0.28094879 -0.00539917 -0.16394293
#> [2,]  0.26060260  0.25201631 -0.13602064  0.04512750 -0.07791046
#> [3,]  0.28094879 -0.13602064  0.38069759 -0.46180782 -0.04561306
#> [4,] -0.00539917  0.04512750 -0.46180782  1.35741376  0.35778699
#> [5,] -0.16394293 -0.07791046 -0.04561306  0.35778699  0.62306807
#> 
#> $Q
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  2.0808185 -0.4417189  0.7283594  1.1373614 -0.5770353
#> [2,] -0.4417189  0.3958207 -0.6479921  0.0219372 -0.1799568
#> [3,]  0.7283594 -0.6479921  3.4900923 -0.6375614 -0.3533824
#> [4,]  1.1373614  0.0219372 -0.6375614  1.6902108 -0.1153048
#> [5,] -0.5770353 -0.1799568 -0.3533824 -0.1153048  1.6097728
#> 

# we could also compute these grammians seperately 
all.equal(gr$P, grammians(obj,'ctr'))
#> [1] TRUE
all.equal(gr$Q, grammians(obj,'obs'))
#> [1] TRUE

# create a random (3 by 3) rational matrix, 
# with a stable and minimum phase statespüace realization
# Note: for the choice "minimum phase" the rational matrix 
# must be square and invertible.
obj = test_stsp(dim = c(3,3), s = 5, bpoles = 1, bzeroes = 1)
gr = grammians(obj, which = 'minimum phase')
gr
#> $P
#>            [,1]        [,2]       [,3]        [,4]       [,5]
#> [1,]  6.8739393  0.60106258 -4.0323640 -2.92750660  1.1841927
#> [2,]  0.6010626  2.00638174 -0.2088486 -0.04179684 -0.2755599
#> [3,] -4.0323640 -0.20884863  3.4931349  2.12353971 -1.2280120
#> [4,] -2.9275066 -0.04179684  2.1235397  2.48924383 -1.1148332
#> [5,]  1.1841927 -0.27555990 -1.2280120 -1.11483316  1.1548715
#> 
#> $Q
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,]  5.490229  4.133327 -5.126688  6.722807  3.207648
#> [2,]  4.133327  4.981322 -4.954169  6.768255  4.164107
#> [3,] -5.126688 -4.954169  8.840520 -8.468811 -4.146512
#> [4,]  6.722807  6.768255 -8.468811 11.933019  4.118414
#> [5,]  3.207648  4.164107 -4.146512  4.118414  6.643613
#> 

# we could also compute these grammians seperately 
all.equal(gr$P, grammians(obj,'ctr'))
#> [1] TRUE
all.equal(gr$Q, grammians(obj,'obs_inv'))
#> [1] TRUE
```
