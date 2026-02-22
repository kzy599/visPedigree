## Test environments
* local macOS Sequoia 15.4, R 4.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  - Version contains large components (1.0.1)
  - Days since last update: ...

This is a major version update from 0.7.1 (last CRAN release) to 1.0.1.

## Changes since last CRAN version (0.7.1)

### New features (v1.0.0)
- `pedmat()`: Comprehensive genetic relationship matrix calculation (A, D, AA, their inverses, and inbreeding coefficients) with Rcpp/Armadillo backend.
- `vismat()`: Relationship matrix visualization with heatmaps and histograms.
- `splitped()`: Detect and split disconnected pedigree components.
- `compact = TRUE` mode for `pedmat()`: Dramatically reduces matrix dimensions for pedigrees with large full-sib families.
- `query_relationship()` and `expand_pedmat()` for querying and restoring compact matrices.
- `tidyped()` now assigns `Family` column automatically; `summary.tidyped()` includes family statistics.
- API standardization: `pedmatrix` → `pedmat`, `n_threads` → `threads`.

### Bug fixes (v1.0.1)
- Fixed compact matrix correctness for parent-offspring pairs.
- Fixed `expand_pedmat()` sibling off-diagonal calculation.
- Fixed `tidyped(..., genmethod = "bottom")` generation alignment.
- Fixed `visped()` edge highlighting and layout with `showf = TRUE`.

## Downstream dependencies
None.
