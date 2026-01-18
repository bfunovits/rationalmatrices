# A Collection of Tools for Rational Matrices

The package `rationalmatrices` collects classes, methods and tools for
creating and manipulating rational matrices, i.e. matrices whose entries
are rational functions.

## Classes

There exist many different representations for rational matrices. This
package in particular deals with "left matrix fraction" representations
(implemented as
[`lmfd`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
class) and "statespace" representations (implemented as
[`stsp`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
class). As a special case of course also polynomial matrices
([`polm`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
class) are implemented. The coefficients of the power series expansion
of the rational function are stored as
[`pseries`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
objects and a collection of values of the rational function may be
stored as
[`zvalues`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
objects.

The package in particular offers tools to convert one representation
into another (equivalent) representation, see e.g.
[`pseries2stsp`](https://bfunovits.github.io/rationalmatrices/reference/pseries2stsp.md)
and
[`pseries2lmfd`](https://bfunovits.github.io/rationalmatrices/reference/pseries2lmfd.md).

## Methods

The package provides "standard" matrix - tools:

- General methods, like [`print`](https://rdrr.io/r/base/print.html),
  [`dim`](https://rdrr.io/r/base/dim.html),
  [`plot`](https://rdrr.io/r/graphics/plot.default.html), ...

- Arithmetic operations, like addition and multiplication. See in
  particular
  [`Ops.ratm`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
  and
  [`%r%`](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md).

- Extract parts of the matrix, transposition, (column or row) bind two
  or more matrices, ...

Some specific methods for rational matrices are

- Compute poles and zeroes of rational matrices.

- Check properties of the rational matrix, like stability and inverse
  stability.

- Some support for the echelon canonical form, see e.g.
  [Kronecker-Indices](https://bfunovits.github.io/rationalmatrices/reference/Kronecker-Indices.md).

- Normal Forms for polynomial matrices, like the Hermite normal form and
  the Smith form. See
  [`hnf`](https://bfunovits.github.io/rationalmatrices/reference/hnf.md),
  [`snf`](https://bfunovits.github.io/rationalmatrices/reference/snf.md),
  [`whf`](https://bfunovits.github.io/rationalmatrices/reference/whf.md)
  and
  [`col_reduce`](https://bfunovits.github.io/rationalmatrices/reference/col_reduce.md).
  Check for "left (co)primeness" with
  [`is.coprime`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md).

- Compute the derivative (with respect to \\z\\).

- Tools related to statespace representations, e.g. computation of
  controllability and observability Grammians and the computation of
  balanced (truncated) realizations. See
  [`grammians`](https://bfunovits.github.io/rationalmatrices/reference/grammians.md),
  [`balance`](https://bfunovits.github.io/rationalmatrices/reference/balance.md),
  ...

## Author(s)

Wolfgang Scherrer and Bernd Funovits

Maintainer: \<bernd.funovits@gmail.com\>

References: See Also: Examples:

## See also

Useful links:

- <https://bfunovits.github.io/rationalmatrices/>

## Author

**Maintainer**: Bernd Funovits <bernd.funovits@gmail.com>
([ORCID](https://orcid.org/0000-0002-8247-6840))

Authors:

- Wolfgang Scherrer <wolfgang.scherrer@tuwien.ac.at>
