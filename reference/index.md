# Package index

## Basic Classes

Different representations of rational matrix functions and special
cases.

- [`lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/lmfd.md)
  : Constructor for Left Matrix Fraction Descriptions (LMFDs)
- [`rmfd()`](https://bfunovits.github.io/rationalmatrices/reference/rmfd.md)
  : Constructor for Right Matrix Fraction Descriptions (RMFDs)
- [`stsp()`](https://bfunovits.github.io/rationalmatrices/reference/stsp.md)
  : Constructor for Statespace Realizations
- [`polm()`](https://bfunovits.github.io/rationalmatrices/reference/polm.md)
  : Constructor for Polynomial Matrices
- [`lpolm()`](https://bfunovits.github.io/rationalmatrices/reference/lpolm.md)
  : Constructor for Laurent Polynomial Matrices
- [`pseries()`](https://bfunovits.github.io/rationalmatrices/reference/pseries.md)
  : Power series Parameters

## Transformation between Different Representations

as.\* methods

- [`as.lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/as.lmfd.md)
  : Coerce to Left Matrix Fraction Description
- [`as.rmfd()`](https://bfunovits.github.io/rationalmatrices/reference/as.rmfd.md)
  : Coerce to Right Matrix Fraction Description
- [`as.stsp()`](https://bfunovits.github.io/rationalmatrices/reference/as.stsp.md)
  : Coerce to Statespace Realization
- [`as.polm()`](https://bfunovits.github.io/rationalmatrices/reference/as.polm.md)
  : Coerece to polynomial object
- [`as.lpolm()`](https://bfunovits.github.io/rationalmatrices/reference/as.lpolm.md)
  : Coerce to Laurent polynom object
- [`pseries2lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2lmfd.md)
  : Construct a LMFD Representation from Impulse Response
- [`pseries2rmfd()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2rmfd.md)
  : Construct an RMFD Representation from Impulse Response
- [`pseries2stsp()`](https://bfunovits.github.io/rationalmatrices/reference/pseries2stsp.md)
  : Ho-Kalman Realization Algorithm
- [`basis2nu()`](https://bfunovits.github.io/rationalmatrices/reference/Kronecker-Indices.md)
  [`nu2basis()`](https://bfunovits.github.io/rationalmatrices/reference/Kronecker-Indices.md)
  [`pseries2nu()`](https://bfunovits.github.io/rationalmatrices/reference/Kronecker-Indices.md)
  : Tools related to Kronecker indices

## Poles and Zeros

- [`poles()`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md)
  [`zeroes()`](https://bfunovits.github.io/rationalmatrices/reference/poles_and_zeroes.md)
  : Poles and Zeroes
- [`zvalue()`](https://bfunovits.github.io/rationalmatrices/reference/zvalue.md)
  : Frequency Response Function
- [`zvalues()`](https://bfunovits.github.io/rationalmatrices/reference/zvalues.md)
  : Frequency Response Function

## Check Object Methods

- [`is.polm()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  [`is.lpolm()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  [`is.lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  [`is.rmfd()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  [`is.stsp()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  [`is.pseries()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  [`is.zvalues()`](https://bfunovits.github.io/rationalmatrices/reference/is.md)
  : Check Objects
- [`is.coprime()`](https://bfunovits.github.io/rationalmatrices/reference/is.coprime.md)
  : Left Prime and Left Coprime Polynomials
- [`is.minimal()`](https://bfunovits.github.io/rationalmatrices/reference/is.minimal.md)
  : Check the Minimality of a Statespace Realization
- [`is.stable()`](https://bfunovits.github.io/rationalmatrices/reference/check.md)
  [`is.miniphase()`](https://bfunovits.github.io/rationalmatrices/reference/check.md)
  : Check Properties of Rational Matrices

## Helpers for Creating Test Objects

- [`test_lmfd()`](https://bfunovits.github.io/rationalmatrices/reference/test_lmfd.md)
  : Create Test Rational Matrix in LMFD Form
- [`test_rmfd()`](https://bfunovits.github.io/rationalmatrices/reference/test_rmfd.md)
  : Create Test Rational Matrix in RMFD Form
- [`test_stsp()`](https://bfunovits.github.io/rationalmatrices/reference/test_stsp.md)
  : Create Test Rational Matrix in Statespace Form
- [`test_polm()`](https://bfunovits.github.io/rationalmatrices/reference/test_polm.md)
  : Create Test Polynomial Matrix
- [`test_lpolm()`](https://bfunovits.github.io/rationalmatrices/reference/test_lpolm.md)
  : Create Test Laurent Polynomial Matrix
- [`test_array()`](https://bfunovits.github.io/rationalmatrices/reference/test_array.md)
  : Create Test Array

## Arithmetic (incl derivative and transpose)

- [`` `%r%` ``](https://bfunovits.github.io/rationalmatrices/reference/ratm_mult.md)
  : Matrix Multiplication of Rational Matrices
- [`Ops(`*`<ratm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/Ops.ratm.md)
  : Arithmetic Ops Group Methods for Rational Matrices
- [`polm_div()`](https://bfunovits.github.io/rationalmatrices/reference/polm_div.md)
  : Division Algorithm for Polynomial Matrices
- [`upgrade_objects()`](https://bfunovits.github.io/rationalmatrices/reference/upgrade_objects.md)
  : Upgrade Objects to Common Class
- [`derivative()`](https://bfunovits.github.io/rationalmatrices/reference/derivative.md)
  : Derivative of a rational Matrix
- [`Ht()`](https://bfunovits.github.io/rationalmatrices/reference/Ht.md)
  : Hermitean Transpose
- [`t(`*`<polm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  [`t(`*`<lpolm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  [`t(`*`<lmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  [`t(`*`<rmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  [`t(`*`<stsp>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  [`t(`*`<pseries>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  [`t(`*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/transpose.md)
  : Rational Matrix Transpose

## Generic Functions

Standard methods for basic S3 classes and their helpers.

- [`` `[`( ``*`<polm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `[`( ``*`<lpolm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `[`( ``*`<lmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `[`( ``*`<rmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `[`( ``*`<stsp>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `[`( ``*`<pseries>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `[`( ``*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `$`( ``*`<lmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `$`( ``*`<rmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `$`( ``*`<stsp>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  [`` `$`( ``*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/extract.md)
  : Extract Parts of a Rational Matrix
- [`` `[<-`( ``*`<polm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/replace.md)
  [`` `[<-`( ``*`<lpolm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/replace.md)
  : Replace Parts of a (Laurent) Polynomial Matrix
- [`rbind(`*`<ratm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/bind.md)
  [`cbind(`*`<ratm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/bind.md)
  : Combine Rational Matrices by Rows or Columns
- [`dim(`*`<polm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  [`dim(`*`<lpolm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  [`dim(`*`<lmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  [`dim(`*`<rmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  [`dim(`*`<stsp>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  [`dim(`*`<pseries>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  [`dim(`*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/dim.md)
  : Dimensions of Objects
- [`str(`*`<lpolm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  [`str(`*`<polm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  [`str(`*`<lmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  [`str(`*`<rmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  [`str(`*`<stsp>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  [`str(`*`<pseries>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  [`str(`*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/str.md)
  : Display the Structure of Objects
- [`print(`*`<lpolm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  [`print(`*`<polm>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  [`print(`*`<lmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  [`print(`*`<rmfd>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  [`print(`*`<stsp>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  [`print(`*`<pseries>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  [`print(`*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/print.md)
  : Print Methods
- [`plot(`*`<pseries>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/plot.md)
  [`plot(`*`<zvalues>`*`)`](https://bfunovits.github.io/rationalmatrices/reference/plot.md)
  : Plot Methods
- [`plot_3D()`](https://bfunovits.github.io/rationalmatrices/reference/plot_3D.md)
  : Plot 3D Arrays
- [`zoom_plot()`](https://bfunovits.github.io/rationalmatrices/reference/zoom_plot.md)
  : Zoom and Scroll
- [`as_tex_matrix()`](https://bfunovits.github.io/rationalmatrices/reference/as_tex_matrix.md)
  : TeX Matrix
- [`as_tex_matrixfilter()`](https://bfunovits.github.io/rationalmatrices/reference/as_tex_matrixfilter.md)
  : TeX Matrix Polynomial Filters
- [`as_tex_matrixpoly()`](https://bfunovits.github.io/rationalmatrices/reference/as_tex_matrixpoly.md)
  : TeX Matrix Polynomials
- [`as_txt_scalarfilter()`](https://bfunovits.github.io/rationalmatrices/reference/as_txt_scalarfilter.md)
  : Coerce Scalar Polynomial Filters to Character Strings
- [`as_txt_scalarpoly()`](https://bfunovits.github.io/rationalmatrices/reference/as_txt_scalarpoly.md)
  : Coerce Scalar Polynomials to Character Strings

## Factorisations, including Wiener-Hopf-Factorisation (WHF) specific functions

- [`whf()`](https://bfunovits.github.io/rationalmatrices/reference/whf.md)
  : Wiener-Hopf Factorization
- [`polm2fwd()`](https://bfunovits.github.io/rationalmatrices/reference/polm2fwd.md)
  : Transforms to Polynomial in Forward Shift
- [`get_fwd()`](https://bfunovits.github.io/rationalmatrices/reference/get_fwd.md)
  [`get_bwd()`](https://bfunovits.github.io/rationalmatrices/reference/get_fwd.md)
  : Forward and Backward Bracket
- [`snf()`](https://bfunovits.github.io/rationalmatrices/reference/snf.md)
  : Smith Normal Form
- [`hnf()`](https://bfunovits.github.io/rationalmatrices/reference/hnf.md)
  : Hermite Normal Form

## Numeric Reduction and Factorisation (including Degree)

- [`col_end_matrix()`](https://bfunovits.github.io/rationalmatrices/reference/col_end_matrix.md)
  : Column End Matrix of a Polynomial Matrix
- [`col_reduce()`](https://bfunovits.github.io/rationalmatrices/reference/col_reduce.md)
  : Construct a Column Reduced Polynomial Matrix
- [`degree()`](https://bfunovits.github.io/rationalmatrices/reference/degree.md)
  : Polynomial Degree
- [`purge_rc()`](https://bfunovits.github.io/rationalmatrices/reference/purge_rc.md)
  : Purge Rows or Columns of a Polynomial Matrix
- [`prune()`](https://bfunovits.github.io/rationalmatrices/reference/prune.md)
  : Prune (Laurent) Matrix Polynomial
- [`schur()`](https://bfunovits.github.io/rationalmatrices/reference/schur.md)
  : Schur Decomposition
- [`ql_decomposition()`](https://bfunovits.github.io/rationalmatrices/reference/ql_decomposition.md)
  [`lq_decomposition()`](https://bfunovits.github.io/rationalmatrices/reference/ql_decomposition.md)
  : QL and LQ Decomposition

## Blaschke et al.

- [`blaschke()`](https://bfunovits.github.io/rationalmatrices/reference/blaschke.md)
  [`blaschke2()`](https://bfunovits.github.io/rationalmatrices/reference/blaschke.md)
  : Blaschke Factors
- [`roots_as_list()`](https://bfunovits.github.io/rationalmatrices/reference/roots_as_list.md)
  : Polynomial Roots as List
- [`reflect_poles()`](https://bfunovits.github.io/rationalmatrices/reference/reflect_poles.md)
  : Reflect Poles of Rational Matrices
- [`reflect_zeroes()`](https://bfunovits.github.io/rationalmatrices/reference/reflect_zeroes.md)
  : Reflect Zeroes of Rational Matrices

## State Space specific Functions

- [`balance()`](https://bfunovits.github.io/rationalmatrices/reference/balance.md)
  : Balanced Realization and Balanced Truncation
- [`companion_matrix()`](https://bfunovits.github.io/rationalmatrices/reference/companion_matrix.md)
  : Companion Matrix of a Polynomial Matrix
- [`ctr_matrix()`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md)
  [`obs_matrix()`](https://bfunovits.github.io/rationalmatrices/reference/ctr_matrix.md)
  : Controllability and Observability Matrix
- [`grammians()`](https://bfunovits.github.io/rationalmatrices/reference/grammians.md)
  : Grammians
- [`state_trafo()`](https://bfunovits.github.io/rationalmatrices/reference/state_trafo.md)
  : State Transformation
- [`lyapunov()`](https://bfunovits.github.io/rationalmatrices/reference/lyapunov.md)
  : Lyapunov Equation

## Block Matrices

- [`bdiag()`](https://bfunovits.github.io/rationalmatrices/reference/bdiag.md)
  : Block Diagonal Matrix
- [`bmatrix()`](https://bfunovits.github.io/rationalmatrices/reference/bmatrix.md)
  : Block matrices
- [`btoeplitz()`](https://bfunovits.github.io/rationalmatrices/reference/btoeplitz.md)
  : Block Toeplitz matrix
- [`bhankel()`](https://bfunovits.github.io/rationalmatrices/reference/bhankel.md)
  : Block Hankel matrix
- [`dbind()`](https://bfunovits.github.io/rationalmatrices/reference/dbind.md)
  : Bind Arrays
- [`array2data.frame()`](https://bfunovits.github.io/rationalmatrices/reference/array2data.frame.md)
  : Coerce arrays to data frames

## Other Stuff

- [`ind2sub()`](https://bfunovits.github.io/rationalmatrices/reference/idx_trafo.md)
  [`sub2ind()`](https://bfunovits.github.io/rationalmatrices/reference/idx_trafo.md)
  : Transform between Linear Index and Matrix Indices
- [`iseq()`](https://bfunovits.github.io/rationalmatrices/reference/iseq.md)
  : Sequence generation
- [`rationalmatrices`](https://bfunovits.github.io/rationalmatrices/reference/rationalmatrices-package.md)
  [`rationalmatrices-package`](https://bfunovits.github.io/rationalmatrices/reference/rationalmatrices-package.md)
  : A Collection of Tools for Rational Matrices
- [`%>%`](https://bfunovits.github.io/rationalmatrices/reference/pipe.md)
  : Pipe operator
