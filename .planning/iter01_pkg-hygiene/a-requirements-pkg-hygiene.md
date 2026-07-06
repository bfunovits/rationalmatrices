---
title: "iter01 — Package Hygiene & Infrastructure (Requirements / Backlog)"
format:
  html:
    embed-resources: true
    page-layout: full
    toc: true
    toc-depth: 3
    toc-expand: 2
---

# iter01 — Package Hygiene & Infrastructure

**Status:** backlog captured — NOT yet started.
**Created:** 2026-07-07 (first iteration of the new workflow, mirroring `acad_RLDM`).
**Source:** `QUALITY_ASSESSMENT.md` (2026-01-18 audit), verified against the current tree at commit `ec4137e`.

This file is the *requirements* artifact for the first iteration of the iteration-folder workflow.
It records the open improvement work identified by the January 2026 quality audit, corrected where the audit is outdated.
Per the workflow rules, no fixes are applied here — research (`b-`), plan (`c-`), execution, changelog (`d-`), and insights (`e-`) follow in this same `iter01_pkg-hygiene/` folder when work begins.

## Corrections to the quality audit (verified 2026-07-07)

The audit predates recent cleanup; three of its findings no longer hold:

- **"Object files committed (8MB)" is outdated** — `git ls-files src/` shows only sources; `src/.gitignore` covers `*.o`/`*.so`. The artifacts (`lyapunov.o`, `RcppExports.o`, `rationalmatrices.so`) exist on disk only. Disk hygiene, not a git problem.
- **"Shiny appears unused" is wrong** — `shiny` is used in `R/09_visualization_tools.R` (interactive 3D visualization). The task is to *document* the dependency, not remove it.
- **"docs/ committed"** — `docs/` is git-ignored; pkgdown deploys via GitHub Actions (`.github/workflows/pkgdown.yaml`).

## Scope and constraints

- **Goal:** bring the package's infrastructure and hygiene to the level of `acad_RLDM` — LICENSE, CI with R CMD check, clean `devtools::check()`, governance docs — then close documentation and test-coverage gaps in priority order.
- **Backward compatibility:** must be preserved — the public API is frozen for this iteration (RLDM and svarmawhf depend on it, including the C++ header `inst/include/rationalmatrices_lyapunov.h`).
- **Out of scope:** new features, new representations, performance optimization, refactoring the large algorithmic functions (`purge_rc`, `col_reduce`, `is.coprime`) beyond documentation.

## Blockers (must fix first)

- [ ] **No LICENSE file** — GPL-3 is declared in `DESCRIPTION` but no `LICENSE`/`LICENSE.md` exists in the root. **DECIDED (2026-07-07):** mirror RLDM exactly — keep `License: GPL-3` in `DESCRIPTION`, add a short `LICENSE.md` pointing to https://www.gnu.org/licenses/gpl-3.0.html, and add `^LICENSE\.md$` to `.Rbuildignore`.
- [ ] **Unknown `devtools::check()` baseline** — no recorded recent full check. Run it, record the result in the changelog, and fix anything above the 0-errors/≤1-warning bar before other work.

## Build / check hygiene

- [ ] Clean compiled artifacts from disk (`src/*.o`, `src/*.so`) and confirm a fresh `devtools::load_all()` rebuild works.
- [ ] Extend `.Rbuildignore` for non-package files so the hidden/extraneous-files NOTE stays clean: `^\.claude$`, `^\.planning$`, `^CLAUDE\.md$`, `^QUALITY_ASSESSMENT\.md$`, `^new-iter\.sh$`, `^z_ignore_this$`, `^.*\.code-workspace$`, `^README\.md$` if it contains non-CRAN badges/links that trigger NOTEs (verify first).
- [ ] Verify `.gitignore` covers `.planning/**/*.html` (rendered Quarto output; the `.md` source is the artifact — RLDM convention).

## CI/CD (highest-impact gap: audit scores CI/CD 3/10)

- [ ] Add `R-CMD-check.yaml` GitHub Actions workflow (r-lib/actions; at minimum Ubuntu + release R, ideally the standard 3-OS matrix). Note: `rationalmatrices` has no GitHub-only dependencies, so the standard template should work directly.
- [ ] Add test-coverage workflow with `covr` + codecov (audit estimates 60–70% coverage; get a real number first).
- [ ] Add a lint configuration (`.lintr`) — low priority, only if it doesn't create noise against the existing style.

## Governance

- [ ] `CONTRIBUTING.md` — mirror RLDM's version, adapted (development setup, numeric-prefix file system, test conventions, Rcpp workflow).
- [ ] README: add badges (R-CMD-check, coverage), a quick-start section, and an installation/development section.
- [ ] ~~`CODE_OF_CONDUCT.md`~~ — **DECIDED (2026-07-07): dropped.** Two-author academic package; revisit only if external contributions appear.

## Documentation

- [ ] Mark the ~18 undocumented internal helpers (from the audit list: `diag_object`, `emult_*`, `expand_letters`, `get_z`, `is_large`, `is_small`, `l2_norm`, `lyapunov_*_cpp`, `nu2stsp_template`, `poly_div`, `poly_rem`, `validate_*`, `whf_scalar`) with `@keywords internal` or `@noRd` as appropriate, then `devtools::document()`.
- [ ] Document the `%r%` operator properly (audit: undocumented despite being the package's central operator).
- [ ] Audit `\dontrun{}` blocks (28 Rd files contain them); promote to runnable examples where feasible.
- [ ] Document the `shiny` dependency (used only by `zoom_plot()` in `R/09_visualization_tools.R`, which already has a `requireNamespace()` guard). **DECIDED (2026-07-07): stays in `Imports`** — the maintainer is the main user and prefers no `DESCRIPTION` change. Just document it.
- [ ] ~~Quick-start vignette~~ — **DECIDED (2026-07-07): deferred.** The README quick-start example (Governance section) covers the discovery need; the two existing vignettes remain the deep documentation.

## Test coverage

Current: 15 test files, ~1,200 lines. Under-tested areas (audit + file inspection):

- [ ] **Visualization** (`R/09_*`) — no dedicated tests at all. At minimum, smoke tests that plot functions run without error on standard objects (use `pdf(NULL)` device).
- [ ] **State-space methods** (`R/06_*`) — grammians, balancing, balanced truncation, controllability/observability need direct tests.
- [ ] **Conversions** — round-trip tests across the representation hierarchy (`polm → stsp → pseries → ...` and back where defined).
- [ ] Edge cases: empty matrices (0 rows/cols), zero polynomial (degree -1), scalar (1×1) cases.

## Deferred tech debt (record dispositions in `e-insights-*.md`)

- [ ] **Large functions** — `purge_rc` (~600 lines), `col_reduce` (~240), `is.coprime` (~210). Do NOT refactor this iteration; add internal documentation (section comments) only if touched anyway.
- [ ] **Downstream contract** — `inst/include/rationalmatrices_lyapunov.h` is included directly by RLDM. Decide whether/how to version this interface (comment header with version + date is probably enough).
- [ ] **CRAN readiness** — the audit estimates 15–20h to CRAN-ready. Not a goal now; note remaining gaps at iteration close.

## Definition of done (for this iteration)

1. `LICENSE` file exists and `devtools::check()` passes clean (0 errors, ≤1 warning), reproducible from a fresh clone.
2. GitHub Actions run R CMD check (and ideally coverage) on push; badges in README.
3. `CONTRIBUTING.md` exists; README has quick-start + development sections.
4. Undocumented-helper NOTE resolved; `%r%` documented.
5. New tests for visualization smoke tests and state-space methods; coverage number recorded.
6. Every remaining backlog item above is either resolved or explicitly deferred with a one-line rationale in `e-insights-*.md`.
