# Controllability and Observability Matrix

The controllability matrix of a statespace realisation
\\k(z)=C(I-Az^{-1})^{-1}B + D\\ is the matrix \$\$
\[B,AB,\dots,A^{o-1}B\] \$\$ and the observability matrix is \$\$
\[C',A'C',\dots,(A')^{o-1}C'\]' \$\$

## Usage

``` r
ctr_matrix(A, B, o = NULL)

obs_matrix(A, C, o = NULL)
```

## Arguments

- A:

  either a
  [`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  object or a square \\(s,s)\\ dimensional matrix.

- B:

  \\(s,n)\\ dimensional matrix. This argument is ignored if `A` is a
  `stsp` object.

- o:

  (non negative) integer. The default value is \\o=s\\.

- C:

  \\(m,s)\\ dimensional matrix. This argument is ignored if `A` is a
  `stsp` object.

## Value

Controllability or observability matrix.
