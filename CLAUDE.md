# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

`rationalmatrices` is an R package providing classes and methods for rational matrices (matrices whose entries are rational functions). It serves as the foundation for the `RLDM` (Rational Linear Dynamic Models) package.

## Build and Development Commands

```bash
# Install dependencies and build
R -e "devtools::install_deps(dependencies = TRUE)"
R -e "devtools::install()"

# Run all tests
R -e "devtools::test()"

# Run a single test file
R -e "testthat::test_file('tests/testthat/test-arithmetic.R')"

# Regenerate documentation (roxygen2)
R -e "roxygen2::roxygenise()"

# Full package check
R -e "devtools::check()"

# Build vignettes
R -e "devtools::build_vignettes()"

# Build pkgdown website
R -e "pkgdown::build_site()"
```

## Core Class Hierarchy

The package implements 7 S3 classes representing different rational matrix representations. Understanding their relationships is essential:

| Class | Description | Internal Storage |
|-------|-------------|------------------|
| `polm` | Polynomial matrices: a(z) = a₀ + a₁z + ... + aₚzᵖ | (m,n,p+1) array |
| `lpolm` | Laurent polynomials with negative powers | (m,n,q+p+1) array + `min_deg` attribute |
| `lmfd` | Left matrix fraction: a⁻¹(z)b(z) | matrix + `order` attribute |
| `rmfd` | Right matrix fraction: d(z)c⁻¹(z) | matrix + `order` attribute |
| `stsp` | State-space: C(z⁻¹I - A)⁻¹B + D | (s+m, s+n) matrix + `order` attribute |
| `pseries` | Power series coefficients | (m,n,lag.max+1) array |
| `zvalues` | Function values at specific z points | (m,n,k) array |

All classes inherit from superclass `ratm`.

**Type coercion rules:**
- Binary operations automatically upgrade to the "higher" representation
- `lpolm` is isolated - cannot be converted to `lmfd`, `rmfd`, `stsp`, or `pseries`
- Use `as.stsp()`, `as.lmfd()`, etc. for explicit conversion

## Key Source Files

- `R/classes.R` - Class constructors (`polm()`, `lmfd()`, `rmfd()`, `stsp()`, etc.)
- `R/arithmetic_methods.R` - Addition, multiplication, `%r%` operator for matrix multiplication
- `R/polm_methods.R` - Normal forms (HNF, SNF, WHF), column reduction, polynomial division
- `R/stsp_methods.R` - Realization algorithms, grammians, controllability/observability
- `R/as_methods.R` - Type coercion between all class types
- `R/tools.R` - Helper functions for matrix operations
- `src/lyapunov.cpp` - C++ Lyapunov solver (RcppArmadillo)

## Important Conventions

- **Matrix multiplication:** Use `%r%` operator for rational matrix multiplication (not `%*%`)
- **Extracting components:** Use `$` accessor (e.g., `lmfd_obj$a`, `lmfd_obj$b`, `stsp_obj$A`)
- **Zero polynomial degree:** The degree of the zero polynomial is -1 (see `degree()`)
- **3D arrays:** Polynomial coefficients stored as (rows, cols, degree+1) arrays

## Testing

Tests are in `tests/testthat/`. Key test files:
- `test-arithmetic.R` - All arithmetic operations
- `test-lyapunov.R` - C++ Lyapunov solver
- `test-coprime.R` - Coprimeness checking
- `test-reflect-poles-zeroes.R` - Pole/zero reflection algorithms

## C++ Code

The package uses Rcpp/RcppArmadillo for performance-critical code:
- `src/lyapunov.cpp` - Lyapunov equation solver
- Auto-generated files: `src/RcppExports.cpp`, `R/RcppExports.R` (do not edit)
- Headers exported via `inst/include/` for dependent packages

## Dependencies

Key imports: `Rcpp`, `RcppArmadillo`, `QZ` (generalized Schur decomposition), `magrittr` (pipe operator)
