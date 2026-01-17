# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Purpose

`rationalmatrices` is an R package for **working with rational matrices** - matrices whose entries are rational functions of z (typically z⁻¹ in time-series/control theory).

### Core Problem

A rational function can be represented in many equivalent forms:
- **Polynomial**: a(z) = a₀ + a₁z + ... + aₚzᵖ
- **Laurent polynomial**: allows negative powers (z⁻¹, z⁻², etc.)
- **Matrix fractions**: a⁻¹(z)b(z) (left) or d(z)c⁻¹(z) (right)
- **State-space**: C(z⁻¹I - A)⁻¹B + D (common in control theory)
- **Impulse response**: h₀ + h₁z⁻¹ + h₂z⁻² + ... (time-domain coefficients)
- **Frequency response**: values evaluated on the unit circle

**The package enables seamless conversion between all these representations**, making it easy to choose the most convenient form for each task. It serves as the foundation for `RLDM` (Rational Linear Dynamic Models).

## Main Functional Areas (Sub-Tasks)

The package is organized around these major task categories:

### 1. **Representation Management** (`R/01_representation_classes.R`, `R/01_representation_conversions.R`)
   - Create rational matrices in different forms: `polm()`, `lpolm()`, `lmfd()`, `rmfd()`, `stsp()`, `pseries()`, `zvalues()`
   - Convert between representations: `as.stsp()`, `as.lmfd()`, `as.rmfd()`, `as.polm()`
   - Validate objects: `is.polm()`, `is.stsp()`, `is.lmfd()`, etc.

### 2. **Realization Algorithms** (`R/01_representation_conversions.R`, `R/02_realization_tools.R`, `vignettes/`)
   - **Hankel-based methods**: Derive LMFD/state-space from impulse response using Kronecker indices
   - **Ho-Kalman algorithm**: State-space realization from power series coefficients
   - **Echelon canonical forms**: Structured representations for analysis

### 3. **Arithmetic Operations** (`R/03_arithmetic_operations.R`)
   - Basic: addition (`+`), subtraction (`-`), element-wise multiplication (`*`)
   - Matrix multiplication: `a %r% b` (special rational operator)
   - Power: `a^k` for both positive and negative integers
   - Binding: `rbind()`, `cbind()` to combine matrices

### 4. **Polynomial Manipulation** (`R/04_polynomial_methods.R`)
   - Normal forms: **HNF** (Hermite), **SNF** (Smith), **WHF** (Wiener-Hopf)
   - Column reduction and echelon forms
   - Polynomial division: `a // b`, `a %% b`
   - Degree computation, column end matrices

### 5. **Analysis & Properties** (`R/05_analysis_poles_zeroes.R`, `R/is_methods.R`)
   - **Poles/Zeros**: Compute roots via eigenvalue decomposition
   - **Stability**: Check if all poles inside unit circle
   - **Minimality**: Verify minimal state-space dimension via Hankel matrix rank
   - **Coprimeness**: Test if polynomial pair has common factors (via singular pencils)

### 6. **State-Space Tools** (`R/06_statespace_methods.R`)
   - **Grammians**: Controllability and observability Grammians
   - **Balancing**: Transform to balanced realization (diagonalize Grammians)
   - **Balanced truncation**: Approximate high-order systems
   - **Controllability/Observability matrices**: System property analysis

### 7. **Numerical Methods** (`src/lyapunov.cpp`, `inst/include/rationalmatrices_lyapunov.h`)
   - **Lyapunov equation solver**: P = APA' + Q (Schur decomposition method)
   - **Stability verification**: Check A eigenvalues
   - C++ header-only implementation for speed and downstream use

### 8. **Pole/Zero Reflection** (`R/08_reflection_poles_zeroes.R`)
   - **Allpass matrices**: Multiply to flip poles/zeros
   - **Blaschke factors**: Univariate and multivariate pole/zero reflection
   - Maintain spectral properties (H2 norm, spectral density)

### 9. **Visualization** (`R/09_visualization_methods.R`, `R/09_visualization_tools.R`)
   - **Frequency response plots**: Magnitude, phase, Nyquist diagrams
   - **Pole-zero diagrams**: Show roots on complex plane
   - **3D visualization**: For matrix-valued functions

### 10. **Utilities & Helpers** (`R/10_utilities_tools.R`, `R/10_utilities_munkres.R`)
   - Matrix operations: Transpose `t()`, Hermitian transpose `Ht()`
   - Derivatives: `derivative()` with respect to z
   - Kronecker products, Hungarian algorithm for matching
   - Array/vector manipulation helpers

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

The package implements 7 S3 classes representing different rational matrix representations. All inherit from `ratm` (rational matrix).

### Representation Classes (Type Hierarchy)

| Class | Best For | Internal Storage | Notes |
|-------|----------|------------------|-------|
| `polm` | Polynomial coefficients a(z) = a₀ + a₁z + ... + aₚzᵖ | (m,n,p+1) array | Foundation; lowest in hierarchy |
| `lpolm` | Laurent polynomials with negative powers z⁻¹, z⁻², ... | (m,n,q+p+1) array + `min_deg` | **Isolated**: can't convert to other MFD/stsp types |
| `lmfd` | **Left** matrix fraction a⁻¹(z)b(z) - for analysis | matrix + `order` | Common in control theory |
| `rmfd` | **Right** matrix fraction d(z)c⁻¹(z) | matrix + `order` | Dual to LMFD; transpose relationship |
| `stsp` | State-space C(z⁻¹I-A)⁻¹B+D - for computation | (s+m, s+n) matrix | Most operations optimized for this |
| `pseries` | Impulse response h₀ + h₁z⁻¹ + h₂z⁻² + ... | (m,n,lag+1) array | Power series coefficients; used for realization |
| `zvalues` | Frequency response: function values at z points | (m,n,k) array | Evaluated on unit circle or arbitrary z |

### Type Coercion System

**Automatic coercion hierarchy** (for binary operations like `a + b`):
```
matrix ≺ polm ≺ lmfd ≺ stsp ≺ pseries ≺ zvalues
```

Binary operations automatically upgrade both operands to their maximum type:
- `polm + lmfd` → result is `lmfd`
- `stsp + pseries` → result is `pseries`

**Important rules:**
- `lpolm` is **isolated** - cannot automatically convert to `lmfd`, `rmfd`, `stsp`, or `pseries`
  - Use explicit `as.stsp(as.polm(lpolm_obj))` if needed
- `lmfd` objects in operations are promoted to `stsp` (control theory preference)
- Use `as.stsp()`, `as.lmfd()`, `as.rmfd()`, etc. for explicit conversion

## Key Source Files

- `R/01_representation_classes.R` - Class constructors (`polm()`, `lmfd()`, `rmfd()`, `stsp()`, etc.)
- `R/03_arithmetic_operations.R` - Addition, multiplication, `%r%` operator for matrix multiplication
- `R/04_polynomial_methods.R` - Normal forms (HNF, SNF, WHF), column reduction, polynomial division
- `R/06_statespace_methods.R` - Realization algorithms, grammians, controllability/observability
- `R/01_representation_conversions.R` - Type coercion between all class types
- `R/10_utilities_tools.R` - Helper functions for matrix operations
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
- `src/lyapunov.cpp` - Rcpp wrapper for Lyapunov equation solver
- `inst/include/rationalmatrices_lyapunov.h` - **Header-only implementation** of the core Lyapunov solver (schema: Schur decomposition-based recursive algorithm). This allows dependent packages like RLDM to directly include and call the C++ solver while preserving pass-by-reference semantics and avoiding SEXP serialization overhead.
- Auto-generated files: `src/RcppExports.cpp`, `R/RcppExports.R` (do not edit)
- Other exported headers: `inst/include/rationalmatrices.h`, `inst/include/rationalmatrices_RcppExports.h`

### Lyapunov Solver Architecture

The Lyapunov solver (`A*P*A' + Q = P`) is implemented via Schur decomposition:
1. Transform A to upper triangular form via Schur decomposition
2. Transform Q to Schur basis
3. Recursively solve for P by iterating from bottom-right to top-left
4. Stability check: all eigenvalues (diagonal of Schur form) must have magnitude < 1
5. Transform solution back to original basis

The core loop is in `solve_lyapunov_core_loop()` in the header file.

## Dependencies

Key imports: `Rcpp`, `RcppArmadillo`, `QZ` (generalized Schur decomposition), `magrittr` (pipe operator)

## Recent Changes (January 2026)

- **Header-only Lyapunov solver** - Extracted core Lyapunov algorithm to `inst/include/rationalmatrices_lyapunov.h` as a C++ header-only library. This enables direct C++ inclusion by dependent packages without Rcpp overhead.
- **Refactored `reflect_poles_zeroes.R`** - Simplified pole/zero reflection logic
- **Simplified `src/lyapunov.cpp`** - Now primarily an Rcpp wrapper around the header-only implementation
- **Added Claude Code configuration** - `.claude/settings.local.json` for local development permissions
