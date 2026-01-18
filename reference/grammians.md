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
#>             [,1]        [,2]        [,3]        [,4]       [,5]
#> [1,]  1.07739784  0.50538423 -0.05728009  0.12505820  0.4495258
#> [2,]  0.50538423  1.87761622 -0.03424442 -0.95850660  1.7004741
#> [3,] -0.05728009 -0.03424442  0.26725447  0.04383366 -0.2638833
#> [4,]  0.12505820 -0.95850660  0.04383366  1.05532416 -1.1579120
#> [5,]  0.44952583  1.70047408 -0.26388326 -1.15791195  2.0176224
#> 
#> $Q
#>           [,1]        [,2]       [,3]        [,4]        [,5]
#> [1,] 2.2335401  0.88589098  0.1979587  1.77875686  0.35699456
#> [2,] 0.8858910  1.31852039 -0.3423208 -0.21070878 -0.04261396
#> [3,] 0.1979587 -0.34232083  1.2636024  0.21796903  0.53815948
#> [4,] 1.7787569 -0.21070878  0.2179690  3.14738379 -0.02437655
#> [5,] 0.3569946 -0.04261396  0.5381595 -0.02437655  0.69695215
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
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,]  0.28565575  0.15967416  0.05261297 -0.03437757  0.09951697
#> [2,]  0.15967416  1.11396776 -0.08427646  0.55940469  0.03650182
#> [3,]  0.05261297 -0.08427646  0.64660514 -0.23693600 -0.27926003
#> [4,] -0.03437757  0.55940469 -0.23693600  1.04710872 -0.10090803
#> [5,]  0.09951697  0.03650182 -0.27926003 -0.10090803  0.80821349
#> 
#> $Q
#>            [,1]       [,2]        [,3]        [,4]        [,5]
#> [1,]  3.6810773 -0.6679773 -0.25048839 -0.13408869 -1.81190242
#> [2,] -0.6679773  2.0605165  0.20156272 -0.19714016  0.45571569
#> [3,] -0.2504884  0.2015627  0.51064851  0.06441107 -0.07452393
#> [4,] -0.1340887 -0.1971402  0.06441107  0.27919111 -0.26344144
#> [5,] -1.8119024  0.4557157 -0.07452393 -0.26344144  2.19753901
#> 

# we could also compute these grammians seperately 
all.equal(gr$P, grammians(obj,'ctr'))
#> [1] TRUE
all.equal(gr$Q, grammians(obj,'obs_inv'))
#> [1] TRUE
```
