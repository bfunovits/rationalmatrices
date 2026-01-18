# A Package for Rational Matrices

`rationalmatrices` collects classes, methods and tools for creating and
manipulating rational matrices, i.e. matrices whose entries are rational
functions. It provides the fundamentals for the R package `RLDM`
(Rational Linear Dynamic Models).

## Installation

`remotes::install_github("bfunovits/rationalmatrices")`

## Classes

There exist many different representations for rational matrices.

This package in particular deals with

- **left matrix fraction** representations (ARMA), implemented as `lmfd`
  class,
- **right matrix fraction** representations, implemented as `rmfd`
  class, and
- **state space** representations, implemented as `stsp` class.

As a special case, polynomial matrices, as `polm` class, and Laurent
polynomial matrices, as `lpolm` class, which allow for negative powers
are implemented as well. The latter class is special in the sense that
even though it is a superset of `polm` object, it cannot be coerced to
`stsp`, `lfmd`, or `rmfd` objects and interaction with other classes
does not make sense in many cases.

The coefficients of the power series expansion of the rational function
are stored as `pseries` objects and a collection of values of the
rational function may be stored as `zvalues` objects.

The package offers tools to convert one representation into another
(equivalent) representation. Some of these methods are quite
sophisticated, e.g. `pseries2stsp` is based on the Ho-Kalman Realization
Algorithm.

## General Matrix Methods

The package provides “standard” generic functions for these objects:

- General methods, like `print`, `dim`, `plot`.
- Arithmetic operations, like addition and (matrix) multiplication.
- Extract parts of the matrix, transposition, (column or row) bind two
  or more matrices.

## Special Methods

The package also provides some specific methods for rational matrices:

- Compute **poles and zeroes** of rational matrices.
- Check properties of the rational matrix, like **stability** and
  inverse stability.
- **Normal Forms** for polynomial matrices, like the Hermite normal
  form, the Smith form, and the Wiener-Hopf factorization.
- Check for **left primeness** with
  [`is.coprime()`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md)
- Compute the **derivative** (with respect to $z$)
- Tools related to state space representations, e.g. computation of
  **controllability and observability** matrices.

## Authors

Wolfgang Scherrer and Bernd Funovits
