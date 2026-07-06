# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Project Overview

`rationalmatrices` is an R package for **working with rational matrices** — matrices whose entries are rational functions of z (typically z⁻¹ in time-series/control theory).
It uses Rcpp/RcppArmadillo for performance-critical computations (Lyapunov solver) and serves as the foundation for the sister package `RLDM` (Rational Linear Dynamic Models), which in turn underpins `svarmawhf`.

This is a **learning project** for the maintainer as well as a working package: explanations of advanced R, Rcpp/C++, and mathematical concepts are valued, especially in iteration insights (see Development Workflow).

**Core problem the package solves:** a rational function can be represented in many equivalent forms —

- **Polynomial**: a(z) = a₀ + a₁z + ... + aₚzᵖ
- **Laurent polynomial**: allows negative powers (z⁻¹, z⁻², etc.)
- **Matrix fractions**: a⁻¹(z)b(z) (left) or d(z)c⁻¹(z) (right)
- **State-space**: C(z⁻¹I - A)⁻¹B + D (common in control theory)
- **Impulse response**: h₀ + h₁z⁻¹ + h₂z⁻² + ... (time-domain coefficients)
- **Frequency response**: values evaluated on the unit circle

**The package enables seamless conversion between all these representations**, making it easy to choose the most convenient form for each task.

# Development Workflow

**Understand requirements → Discuss → Research → Plan → Execute → Verify → Log → Insights**

Development flow is based on subfolders of the form `iter<two-digit-number>_<datestamp|description>` in `.planning/`, e.g. `.planning/iter01_pkg-hygiene` (the datestamp and/or description may be missing).
Use `./new-iter.sh <number> [description]` to create a new iteration folder.
When a Claude Code session starts work on an iteration, ensure the folder is named with its iteration number plus a short description.

All artefacts for an iteration live in the same `iter` subfolder as its requirements file.
If you do research, it goes in `.planning/iter*/b-research-*.md`; the plan goes in `.planning/iter*/c-plan-*.md`; and so on.
(If you start in plan mode, the default plan location may differ — move it into the iter folder.)

| File | Created by | Purpose |
|---|---|---|
| `a-requirements-<desc>.md` | User | Scope and constraints |
| `b-research-<desc>.md` | Agent (if needed) | Ecosystem research, trade-off analysis |
| `c-plan-<desc>.md` | Agent | Atomic task list with verification steps |
| `d-changelog-<desc>.md` | Agent | What changed, starting git commit hash |
| `e-insights-<desc>.md` | Agent | Reflections, learnings, future ideas — IMPORTANT (this is a learning project): include explanations of advanced programming/mathematical concepts and CLAUDE.md suggestions |

## Stages

**Discuss** — Agent reads the requirements file.
If ambiguities exist, ask clarifying questions in the usual Claude Code multiple-choice style (convenient to navigate and submit).
Suggest the best answer and argue why. Resolve before planning.

**Research** (optional) — When the plan needs ecosystem knowledge (R packages, Rcpp/Armadillo patterns, algorithms), research first and save as `b-research-*.md`. Update it if findings change the plan.

**Plan** — Break work into atomic, independently verifiable phases.
Each phase states: files touched, what changes, how to verify. The agent chooses the execution order. Save as `c-plan-*.md`.

**Execute** — Implement one phase at a time.
Run code / tests after each phase to verify before moving on. Do not batch. If a phase fails verification, fix before proceeding.

**Verify / Update docs** — At the end of each iteration, confirm that new or changed public functions/classes are reflected in the Roxygen docs, vignettes, and (if user-facing) the pkgdown site. Run the validation suite (see Task Completion).

**Log** — After all phases complete, write `d-changelog-*.md` (include the starting commit hash).

**Insights** — Write `e-insights-*.md`. Standardized structure for cross-iteration synthesis:

1. Technical learnings (explain the advanced R / Rcpp / mathematical concepts involved — this is a learning project)
2. Process learnings
3. CLAUDE.md suggestions — **as a disposition table** (suggestion | recommended apply/defer/drop | user decision)
4. New skill proposals (with enough detail to create the skill in `.claude/skills/`)
5. Open follow-ups

## Iteration close ritual

Open the insights file with the user, walk the disposition table, and apply accepted CLAUDE.md changes **in the same session** — permission is explicit and fresh.
Append the iteration's technical and process learnings to `.planning/LEARNINGS.md` (the running cross-iteration log) at every close.
Accepted skill proposals are created in the same or next iteration, never deferred beyond that.

## Periodic synthesis

Every ~5 iterations, audit the disposition tables of the last 5 insights files — confirm accepted items landed, re-decide deferred items, reconcile into `.planning/LEARNINGS.md`.

## Multi-phase iterations & model routing

Use sequential phases producing self-contained files, then a synthesis phase that loads all phase files.

| Stage | Model | Rationale |
|---|---|---|
| Algorithm/math design, realization/normal-form logic, hard debugging, synthesis judgment | Capable model (Opus, the default) | Highest-stakes reasoning |
| Routine implementation, refactors, writing tests | Sonnet | Fast, accurate code gen |
| Mechanical work (rendering loops, file moves, doc formatting, categorization) | Haiku | Cost-efficient |

After planning and before executing phases, remind the user to clear context if useful.
Interrupt and flag whenever switching to a cheaper model is possible — **the user is cost sensitive** (especially flag when Haiku would suffice).
State the recommendation/result at phase boundaries.

## Quarto rendering of iteration files

For all agent-written `.md` files in the table above (and `.planning/LEARNINGS.md`), use this YAML header:

```
---
title: <TITLE_STRING>
format:
  html:
    embed-resources: true
    page-layout: full
    toc: true
    toc-depth: 3
    toc-expand: 2
---
```

and render with `quarto render <path_to_md_file>`.
Quarto binary: `/usr/lib/rstudio/resources/app/bin/quarto/bin/quarto` (or the Positron build at `/usr/share/positron/resources/app/quarto/bin/quarto`).
Rendered HTML stays inside the iter folder (it is git-ignored); the `.md` source is the tracked artifact.

## Rules

- **Agent does NOT commit** — the user reviews and commits.
- **Scratch / experimental scripts go in the `iter` folder, prefixed `z_`** — never in the package root, `R/`, or `tests/`. Genuine package unit tests belong in `tests/testthat/` (named after the functionality they test). Tools that prove durable are promoted to `R/` or `inst/` with explicit user approval.
- **DO NOT ADD TO CLAUDE.md WITHOUT EXPLICIT PERMISSION** — the disposition table at iteration close is where permission is requested.
- In each code file (R or C++), include a header comment / Roxygen block stating: last-updated date, a short summary, and the other files it connects to.
- After modifying C++: run `Rcpp::compileAttributes()` then `devtools::load_all()`.
- Ask clarifying questions; do not silently guess about intent, architecture, or requirements.
- Simplest solution first; implement the minimum that works. No unrequested abstractions.
- Flag uncertainty explicitly.

## Handoff files

When work spans sessions, write `.planning/iter*/handoff-<desc>.md` (use the `handoff-writer` skill if available): current context, completed files, remaining tasks, key decisions, environment details.

# Build Commands

```r
# Development workflow
devtools::load_all()           # Load package in dev mode
devtools::document()           # Generate Roxygen docs + NAMESPACE
devtools::test()               # Run testthat tests
devtools::check()              # Full package check (target: 0 errors, <=1 warning)
devtools::build()              # Build source package

# Run a single test file
testthat::test_file("tests/testthat/test-lyapunov.R")

# Rcpp workflow (after modifying C++ code)
Rcpp::compileAttributes()      # Regenerate RcppExports.cpp/R
devtools::load_all()           # Reload with new compiled code

# Documentation site
pkgdown::build_site()          # Build pkgdown site locally (output -> docs/)
devtools::build_vignettes()
```

# Architecture

## Core class hierarchy (S3)

The package implements 7 S3 classes representing different rational matrix representations. All inherit from `ratm` (rational matrix).

| Class | Best For | Internal Storage | Notes |
|-------|----------|------------------|-------|
| `polm` | Polynomial coefficients a(z) = a₀ + a₁z + ... + aₚzᵖ | (m,n,p+1) array | Foundation; lowest in hierarchy |
| `lpolm` | Laurent polynomials with negative powers z⁻¹, z⁻², ... | (m,n,q+p+1) array + `min_deg` | **Isolated**: can't convert to other MFD/stsp types |
| `lmfd` | **Left** matrix fraction a⁻¹(z)b(z) - for analysis | matrix + `order` | Common in control theory |
| `rmfd` | **Right** matrix fraction d(z)c⁻¹(z) | matrix + `order` | Dual to LMFD; transpose relationship |
| `stsp` | State-space C(z⁻¹I-A)⁻¹B+D - for computation | (s+m, s+n) matrix | Most operations optimized for this |
| `pseries` | Impulse response h₀ + h₁z⁻¹ + h₂z⁻² + ... | (m,n,lag+1) array | Power series coefficients; used for realization |
| `zvalues` | Frequency response: function values at z points | (m,n,k) array | Evaluated on unit circle or arbitrary z |

## Type coercion system

**Automatic coercion hierarchy** (for binary operations like `a + b`):

```
matrix ≺ polm ≺ lmfd ≺ stsp ≺ pseries ≺ zvalues
```

Binary operations automatically upgrade both operands to their maximum type:

- `polm + lmfd` → result is `lmfd`
- `stsp + pseries` → result is `pseries`

**Important rules:**

- `lpolm` is **isolated** — cannot automatically convert to `lmfd`, `rmfd`, `stsp`, or `pseries`. Use explicit `as.stsp(as.polm(lpolm_obj))` if needed.
- `lmfd` objects in operations are promoted to `stsp` (control theory preference).
- Use `as.stsp()`, `as.lmfd()`, `as.rmfd()`, etc. for explicit conversion.

## Key file organization (numeric functional system)

R files use a numeric prefix system organized by **purpose/workflow**, not by object type (S3 dispatch handles representation differences).

| Prefix | Category | Content |
|--------|----------|---------|
| **01** | Representations | Class constructors (`polm()`, `lpolm()`, `lmfd()`, `rmfd()`, `stsp()`, `pseries()`, `zvalues()`) and conversions (`as.stsp()`, `as.lmfd()`, ...) |
| **02** | Realization | Hankel-based methods, Ho-Kalman algorithm, echelon canonical forms |
| **03** | Arithmetic | `+`, `-`, `*`, `%r%` (rational matrix multiplication), `^`, `rbind()`/`cbind()`, transpose `t()`/`Ht()` |
| **04** | Polynomial methods | Normal forms (HNF, SNF, WHF), column reduction, polynomial division (`//`, `%%`), `degree()` |
| **05** | Analysis | Poles/zeros, stability, minimality, coprimeness, `derivative()` |
| **06** | State-space tools | Grammians, balancing, balanced truncation, controllability/observability |
| **07** | Numerical | Lyapunov solver (R wrapper + Rcpp), Schur decomposition helpers |
| **08** | Reflection | Allpass matrices, Blaschke factors for pole/zero reflection |
| **09** | Visualization | Frequency response plots, pole-zero diagrams, 3D visualization |
| **10** | Utilities | `dim`, extraction (`$`), `print`, `str`, validation, Munkres/Hungarian algorithm, helpers |

## C++ code (`src/`)

- `src/lyapunov.cpp` — Rcpp wrapper for the Lyapunov equation solver.
- `inst/include/rationalmatrices_lyapunov.h` — **Header-only implementation** of the core Lyapunov solver (Schur decomposition-based recursive algorithm). Dependent packages like RLDM include and call the C++ solver directly, preserving pass-by-reference semantics and avoiding SEXP serialization overhead.
- Other exported headers: `inst/include/rationalmatrices.h`, `inst/include/rationalmatrices_RcppExports.h`.
- Auto-generated files: `src/RcppExports.cpp`, `R/RcppExports.R` (do not edit).

**Lyapunov solver architecture** (`A*P*A' + Q = P`, solved via Schur decomposition):

1. Transform A to upper triangular form via Schur decomposition
2. Transform Q to Schur basis
3. Recursively solve for P by iterating from bottom-right to top-left
4. Stability check: all eigenvalues (diagonal of Schur form) must have magnitude < 1
5. Transform solution back to original basis

The core loop is in `solve_lyapunov_core_loop()` in the header file.

# Key Entry Points

- Construction: `polm()`, `lpolm()`, `lmfd()`, `rmfd()`, `stsp()`, `pseries()`, `zvalues()`
- Conversion: `as.polm()`, `as.lmfd()`, `as.rmfd()`, `as.stsp()`, `as.pseries()`, `as.zvalues()`
- Arithmetic: `%r%` (matrix multiplication), `+`, `-`, `*`, `^`, `t()`, `Ht()`
- Analysis: `poles()`, `zeroes()`, `is.stable()`, `is.minimal()`, `is.coprime()`
- Polynomial: `hnf()`, `snf()`, `whf()`, `col_reduce()`, `degree()`
- State space: `grammians()`, `balance()`, `ctr_matrix()`, `obs_matrix()`
- Reflection: `reflect_poles()`, `reflect_zeroes()`, `blaschke()`

# Important Conventions

- **Matrix multiplication:** Use `%r%` operator for rational matrix multiplication (not `%*%`)
- **Extracting components:** Use `$` accessor (e.g., `lmfd_obj$a`, `lmfd_obj$b`, `stsp_obj$A`)
- **Zero polynomial degree:** The degree of the zero polynomial is -1 (see `degree()`)
- **3D arrays:** Polynomial coefficients stored as (rows, cols, degree+1) arrays

# Known Issues & Fragile Areas

Open work is tracked in `.planning/iter01_pkg-hygiene/a-requirements-pkg-hygiene.md`; the underlying audit is `QUALITY_ASSESSMENT.md` (January 2026, partially outdated — see the requirements file for corrections). Highlights:

- **No LICENSE file** — GPL-3 is declared in `DESCRIPTION` but the file is missing.
- **No R CMD check CI** — `.github/workflows/` only builds pkgdown; no automated check/test/coverage.
- **Large functions:** `purge_rc` (~600 lines), `col_reduce` (~240), `is.coprime` (~210) — mathematically justified but fragile; add tests before changing.
- **Under-tested areas:** visualization (no dedicated tests) and state-space methods.
- **`lpolm` isolation** is by design, not a bug — do not "fix" it.

# Testing

Tests are in `tests/testthat/` (15 files), named after functionality (e.g. `test-lyapunov.R`, `test-coprime.R`, `test-reflect-poles-zeroes.R`).

**Test patterns:**

- Use `testthat::test_file("tests/testthat/test-*.R")` for specific files; `devtools::test(filter = "pattern")` for patterns.
- Common assertions: `expect_equal()`, `expect_true()`, `expect_error()`, `expect_no_error()`.
- Set seeds with `set.seed()` for reproducible stochastic tests (random test matrices via `test_polm()`, `test_lmfd()`, `test_stsp()`, ...).
- Numerical comparisons need tolerances — realization/normal-form algorithms are exact only up to floating point.

# Documentation

## Vignettes

1. `vignettes/a_rational_matrices.Rmd` — main introduction to rational matrices and the class system.
2. `vignettes/b_technical_details_ratm.Rmd` — technical details of the `ratm` classes.

**Vignette file naming:** vignette filenames must start with a letter — `R CMD check` enforces `^[A-Za-z][A-Za-z0-9._-]+$`, so a leading digit triggers a WARNING. Control reading order via `_pkgdown.yml` (`articles:`) and `\VignetteIndexEntry`, not via numeric filename prefixes.

## Documentation build

- `docs/` is the **pkgdown build output** (deployed to https://bfunovits.github.io/rationalmatrices/ via GitHub Actions) — it is generated, not hand-edited, and git-ignored.
- Roxygen2 with markdown support; references managed via Rdpack; `_PACKAGE` convention for package-level docs.

# Task Completion

Run this validation sequence before considering a task done (and at every iteration close):

```r
Rcpp::compileAttributes()      # only if C++ was modified
devtools::load_all()           # package loads cleanly
devtools::document()           # docs + NAMESPACE up to date
devtools::test()               # tests pass
devtools::check()              # 0 errors, <=1 warning
```

Checklist: new R files placed by numeric prefix; exported functions have Roxygen with `@param`/`@return`/`@examples`; internal helpers use `@noRd` or `@keywords internal`; `devtools::run_examples()` passes; new dependencies added to `DESCRIPTION`.

When backgrounding a long command (e.g. `devtools::check()`), wait on its PID with `until ! kill -0 <PID>; do sleep 5; done` — do **not** use `pgrep -f <scriptname>`, which self-matches the wait command's own command line and so never exits.

# Coding Conventions

## Markdown (for `.md`/`.qmd` that Quarto renders)

Quarto requires these blank lines or lists/tables do not render:

- BLANK LINE after every heading (`#`, `##`, or bold text used as an informal header / followed by a colon).
- BLANK LINE before every bullet list, numbered list, and table — even when the preceding line ends with a colon.
- BLANK LINE between consecutive sections.
- Start every sentence on a new line (easier to track in git).
- No harm in extra blank lines — better too many than too few.

## R

**Package code (`R/`, `src/`)** is base R + S3 + Rcpp/Armadillo and existing files use `<-` in places. When editing existing package code, **match the surrounding file's style** for consistency.

**New analysis / iteration scripts and new standalone code** — house style:

- Use ` = ` for assignment (not ` <- `).
- Prefer tidyverse, nested tibbles, functional programming with `purrr` (`map` over loops).
- Use the `logger` package for logging.
- **Never use the `data.table` package.**

**Always (both):** snake_case; meaningful mathematical names; S3 dispatch for the `ratm` classes; Roxygen2 with markdown for exported functions (`@inheritParams` for shared params, Rdpack for references); place new R files per the numeric-prefix system; `@noRd` for internal helpers.

## C++ (`src/`, `inst/include/`)

- C++11+; explicit `arma::` prefix (no `using namespace arma;`).
- Google C++ style: 2-space indentation, ~80-char line limit, snake_case.
- Doxygen-style comments (`/** ... */`) with `@brief`, `@param`, `@return`.
- After any change: `Rcpp::compileAttributes()` then `devtools::load_all()`; never hand-edit `RcppExports.*`.
- Changes to `inst/include/rationalmatrices_lyapunov.h` affect **downstream packages** (RLDM) — keep the header self-contained and re-run downstream checks when the signature changes.

## Quarto

- Don't use a bare `$` (it needs escaping in Quarto) — write USD or escape it.
- For package docs, R code blocks (`{r}`) work because Quarto uses the system R installation.

**Pre-render checklist:** blank lines after headings / before lists; `library()` calls present in each chunk using non-base functions; no bare `$`; all referenced images exist on disk.
