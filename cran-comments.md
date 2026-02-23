## Test environments
* local macOS Sequoia 15.4, R 4.5.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

This is a major version update from 0.7.1 (last CRAN release) to 1.0.1.

## Changes since last CRAN version (0.7.1)

### New features (v1.0.1)
- `pedmat()`: Comprehensive genetic relationship matrix calculation (A, D, AA, their inverses, and inbreeding coefficients) with Rcpp/Armadillo backend.
- `inbreed()`: Re-implemented inbreeding coefficient calculation natively in C++ using the Meuwissen & Luo (1992) algorithm, removing the hard dependency on the `nadiv` package (moved to `Suggests`).
- `vismat()`: Relationship matrix visualization with heatmaps and histograms.
- `splitped()`: Detect and split disconnected pedigree components.
- `compact = TRUE` mode for `pedmat()`: Dramatically reduces matrix dimensions for pedigrees with large full-sib families.
- `query_relationship()` and `expand_pedmat()` for querying and restoring compact matrices.
- `tidyped()` now assigns `Family` column automatically; `summary.tidyped()` includes family statistics.
- `tidyped()`: Added `genmethod` parameter for flexible generation assignment (`"top"` or `"bottom"`).
- API standardization: `pedmatrix` → `pedmat`, `n_threads` → `threads`.

### Bug fixes (v1.0.1)
- Fixed `visped()` edge highlighting and layout with `showf = TRUE`.
- Fixed `tidyped(..., genmethod = "bottom")` generation alignment to prioritize sibling consistency.
- Fixed compact matrix correctness for parent-offspring pairs and sibling off-diagonal calculation in `expand_pedmat()`.

## Downstream dependencies
None.
