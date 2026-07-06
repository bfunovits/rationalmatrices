---
title: "iter01 — Package Hygiene & Infrastructure (Plan)"
format:
  html:
    embed-resources: true
    page-layout: full
    toc: true
    toc-depth: 3
    toc-expand: 2
---

# iter01 — Package Hygiene & Infrastructure (Plan)

**Baseline commit:** `ec4137e` (`dont gitignore Rd files`).

**Requirements:** `a-requirements-pkg-hygiene.md`.

**Goal:** LICENSE + clean `devtools::check()` + R CMD check CI + governance docs, then documentation and test-coverage gaps.
Public API is frozen; no new features; no refactoring of the large algorithmic functions.

This plan breaks the work into atomic, independently verifiable phases.
Blockers first (Phases 1–2), then CI (Phase 3), governance (Phase 4), documentation (Phase 5), tests (Phase 6), tech-debt dispositions + final validation and logs (Phase 7).
Phases 3–6 are independent of each other and can be reordered or split across sessions; each ends with a verification step.

## Model routing (cost-sensitive — flag at each boundary)

| Phase | Suggested model | Rationale |
|---|---|---|
| 1 — LICENSE + disk cleanup | **Haiku** | Pure mechanical file work |
| 2 — baseline check + `.Rbuildignore` | **Haiku** (run check), Sonnet to interpret output | Mechanical edits; reading the check log is routine |
| 3 — CI workflows | Sonnet | Template adaptation, bounded |
| 4 — governance docs | Sonnet | Writing from RLDM templates |
| 5 — documentation | Sonnet | Judgment on `\dontrun` promotion and `shiny` placement, but bounded |
| 6 — tests | Sonnet | Routine test writing against existing patterns |
| 7 — dispositions + validation + logs | **Opus** | Synthesis and tech-debt judgment is high-value |

Interrupt and recommend the cheaper model at the start of Phases 1 and 2.

---

## Phase 1 (BLOCKER) — LICENSE file + `src/` disk cleanup

**Files touched:** `LICENSE.md` (new), `.Rbuildignore`; disk-only deletions in `src/`. `DESCRIPTION` is NOT touched.

**Changes:**

1. Add the license the way RLDM does it (**decided 2026-07-07**): keep `License: GPL-3` in `DESCRIPTION`; create a short `LICENSE.md` containing only a pointer ("This package is licensed under the GPL-3 license. See https://www.gnu.org/licenses/gpl-3.0.html for details."); add `^LICENSE\.md$` to `.Rbuildignore`. Do NOT use `usethis::use_gpl3_license()` (it would change `DESCRIPTION` to `GPL (>= 3)` and diverge from the sibling repos).
2. Delete compiled artifacts from disk: `src/lyapunov.o`, `src/RcppExports.o`, `src/rationalmatrices.so` (already git-ignored via `src/.gitignore` — no ignore-rule change needed).

**Verify:**

- `ls LICENSE.md` exists; `grep License DESCRIPTION` is consistent with it.
- `ls src/` shows only: `.gitignore` (if present), `Makevars`, `Makevars.win`, `lyapunov.cpp`, `RcppExports.cpp`.
- `devtools::load_all()` recompiles cleanly from scratch.
- `git status --porcelain` shows only the intended additions.

---

## Phase 2 (BLOCKER) — Baseline `devtools::check()` + `.Rbuildignore` hygiene

**Files touched:** `.Rbuildignore`; whatever small fixes the check surfaces.

**Changes:**

1. Extend `.Rbuildignore` with: `^\.claude$`, `^\.planning$`, `^CLAUDE\.md$`, `^QUALITY_ASSESSMENT\.md$`, `^new-iter\.sh$`, `^z_ignore_this$`, `^.*\.code-workspace$`, `^LICENSE\.md$` (if the full-text route was taken in Phase 1).
2. Run `devtools::check()` in the background (wait on the PID, not `pgrep`); save the log to `.planning/iter01_pkg-hygiene/z_check_baseline.log`.
3. Triage: fix everything mechanical (NOTEs about files, stale Rd cross-references, example failures). Anything algorithmic goes on the Phase-7 disposition list instead of being hot-fixed.

**Verify:**

- `devtools::check()` result: 0 errors, ≤1 warning; log saved in the iter folder.
- The hidden/extraneous-files NOTE is gone.

---

## Phase 3 — CI: R CMD check + coverage

**Files touched:** `.github/workflows/R-CMD-check.yaml` (new), `.github/workflows/test-coverage.yaml` (new), `README.md` (badges only), `codecov.yml` (optional).

**Changes:**

1. Add the standard `r-lib/actions` R-CMD-check workflow (v2 examples). Start with `ubuntu-latest` + release R; extend to the 3-OS matrix only if the single-OS run is green. System deps: none unusual (Rcpp/RcppArmadillo/QZ compile from source on the runner; setup-pandoc needed for vignettes).
2. Add the `covr`/codecov test-coverage workflow. Record the first real coverage number in the changelog.
3. Add both badges to `README.md`.

**Verify:**

- Workflows are valid YAML (`gh workflow list` after the user pushes, or a dry parse).
- Since the agent does not commit/push: leave verification-on-CI as an explicit checklist item for the user in the changelog; locally run `covr::package_coverage()` once to confirm it executes and to get the number.

---

## Phase 4 — Governance docs

**Files touched:** `CONTRIBUTING.md` (new), `README.md`.

**Changes:**

1. Adapt RLDM's `CONTRIBUTING.md`: development setup (`devtools`, `remotes`), numeric-prefix file placement, test conventions, Rcpp workflow (`compileAttributes()` + `load_all()`), the downstream-header caveat (`inst/include/`), validation sequence.
2. README: add quick-start (create a `polm`, convert to `stsp`, compute poles — a 10-line example), installation (`remotes::install_github("bfunovits/rationalmatrices")`), and a development section pointing to CONTRIBUTING.md. The README quick-start also covers the deferred quick-start vignette (decided 2026-07-07).
3. `CODE_OF_CONDUCT.md`: **decided 2026-07-07 — dropped.** Record the disposition in `e-insights-*.md`; no file created.

**Verify:**

- `devtools::check()` still clean (README/CONTRIBUTING are build-ignored or harmless).
- Quick-start example actually runs (`reprex` or paste into R).

---

## Phase 5 — Documentation gaps

**Files touched:** roxygen blocks across `R/`, regenerated `man/`. `DESCRIPTION` is NOT touched.

**Changes:**

1. Add `@keywords internal` (exported-but-internal) or `@noRd` (unexported) to the ~18 helpers listed in the requirements; regenerate docs.
2. Write proper roxygen for `%r%` (description, `@param`, `@return`, `@examples` showing `polm %r% polm` and mixed-class multiplication with the coercion rule).
3. Audit `\dontrun{}` in the 28 affected Rd-producing blocks: promote to plain `@examples` where they run in <5s and have no side effects; downgrade to `\donttest{}` where slow; keep `\dontrun{}` only for interactive/shiny examples — record counts.
4. `shiny` dependency: **decided 2026-07-07 — stays in `Imports`** (maintainer is the main user; no `DESCRIPTION` change wanted). Only add a note in the `zoom_plot()` roxygen that shiny powers the interactive display; the existing `requireNamespace()` guard at `R/09_visualization_tools.R:910` stays as is.

**Verify:**

- `devtools::document()` clean; undocumented-object NOTE gone.
- `devtools::run_examples()` passes.
- `?"%r%"` renders with examples.
- `git diff DESCRIPTION` is empty (no dependency changes in this phase).

---

## Phase 6 — Test coverage

**Files touched:** new files in `tests/testthat/`: `test-visualization.R`, `test-statespace-methods.R`, `test-conversion-roundtrip.R`, `test-edge-cases.R`.

**Changes:**

1. **Visualization smoke tests** — wrap plot calls in `pdf(NULL)`/`dev.off()`, assert `expect_no_error()` for the main plot methods on small `polm`/`stsp`/`zvalues` objects. Skip shiny-dependent functions with `skip_if_not_installed("shiny")`.
2. **State-space methods** — `grammians()` (solves the Lyapunov equation: verify `P - A P A' - Q ≈ 0`), `balance()` (Grammians diagonal and equal after balancing), balanced truncation (impulse-response error small for a fast-decaying system), `ctr_matrix()`/`obs_matrix()` ranks for known-minimal systems.
3. **Round-trip conversions** — `polm → stsp → pseries` vs direct `pseries(polm)`; `lmfd ↔ stsp`; compare via `pseries` coefficients with tolerance `1e-8`. Use `set.seed()` + the `test_*` generators.
4. **Edge cases** — zero polynomial (degree -1), 1×1 scalars, 0-row/0-column matrices through construction, arithmetic, and printing.

**Verify:**

- `devtools::test()` all green; no test takes >5s.
- `covr::package_coverage()` re-run; record the delta vs the Phase-3 baseline in the changelog.

---

## Phase 7 — Tech-debt dispositions, final validation, logs

**Files touched:** `.planning/iter01_pkg-hygiene/d-changelog-pkg-hygiene.md`, `e-insights-pkg-hygiene.md`, `.planning/LEARNINGS.md`.

**Changes:**

1. Walk the deferred list from the requirements (large functions, header versioning, CRAN readiness, quick-start vignette, `.lintr`) and record a one-line disposition each (fix now / defer with reason / drop).
2. Run the full validation sequence: `Rcpp::compileAttributes()` (if C++ touched — should not be), `devtools::load_all()`, `document()`, `test()`, `check()`.
3. Write `d-changelog-pkg-hygiene.md` (baseline commit `ec4137e`, per-phase summary, check/coverage numbers, the user's push-and-watch-CI checklist).
4. Write `e-insights-pkg-hygiene.md` with the standardized structure, including the CLAUDE.md disposition table; render all iteration `.md` files with Quarto; start `.planning/LEARNINGS.md`.

**Verify:**

- Validation sequence clean (0 errors, ≤1 warning).
- All five artifact files exist in the iter folder; HTML rendered; `.md` tracked.
- Iteration close ritual scheduled with the user.
