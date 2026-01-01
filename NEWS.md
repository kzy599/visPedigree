# Changes in version 0.6.2 released on 01 Jan 2026
## New features
1. Added `summary()` method for `tidyped` objects to provide quick pedigree statistics (number of individuals, founders, sex distribution, etc.).

## Bug fixes
1. Fixed an issue where `tidyped(..., inbreed=TRUE)` failed due to incorrect class assignment order.
2. Fixed `visped(..., showf=TRUE)` to gracefully handle missing `f` columns by warning the user instead of erroring.
3. Fixed broken internal navigation links in package vignettes.

# Changes in version 0.6.1 released on 30 Dec 2025
## New features
1. Implemented opaque highlighting effects for better visualization clarity.
2. Added `trace` option to `visped()` to control ancestry tracing direction.

# Changes in version 0.6.0 released on 28 Dec 2025
## New features
1. Implemented strict S3 class structure for `tidyped` objects with `new_tidyped()` constructor and `validate_tidyped()` validator to ensure data integrity.

# Changes in version 0.5.0 released on 26 Dec 2025
## New features
1. Added `highlight` parameter to `visped()` function. Users can now highlight specific individuals using a character vector of IDs or a list for custom colors.
2. Added `showf` parameter to `visped()` function to display inbreeding coefficients on the pedigree graph.
3. Added `inbreed` parameter to `tidyped()` function to calculate inbreeding coefficients using the `nadiv` package.
4. Refactored `inbreed()` function as a standalone tool that operates on `tidyped` objects.
5. Optimized `repeloverlap()` function using `data.table` for significantly better performance.

## Bug fixes
1. Fixed a critical crash in `visped()` when combining `compact = TRUE`, `highlight`, and `showf = TRUE` by refactoring `ped2igraph()` to delay label modification until after layout calculation.
2. Fixed documentation grammar and phrasing across all functions for CRAN compliance.
3. Fixed `R CMD check` notes related to `data.table` non-standard evaluation by adding `R/globals.R`.

# Changes in version 0.4.1 released on 25 Dec 2025

# Changes in version 0.2.6 released on 31 Mar 2020
## New features
## Bug fixes
1. Fixed a bug that the number of generations for candidates would be traced to n+1 when tracegen=n. This bug is found by Mianyu Liu.

# Changes in version 0.2.5 released on 25 Feb 2020
## New features
## Bug fixes
1. The tidyped() does not work with trace='all' in [certain cases](https://github.com/luansheng/visPedigree/issues/2#issue-568599008)

# Changes in version 0.2.4.1 released on 24 Feb 2020
## New features
## Bug fixes
1. An unexpected column with the name as NA occured when a tidyped object is tidyed again using the tidyped()

# Changes in version 0.2.4 released on 12 June 2019
## New features
## Bug fixes
1. The data.table used as the input parameter 'ped' may be changed in tidyped() and visped().


# Changes in version 0.2.3 released on 05 Mar 2019
## New features
## Bug fixes
1. The generation number of individuals is not inferred rightly.

# Changes in version 0.2.2 released on 28 Jan 2019
## New features
## Bug fixes
1. The tidied pedigree will not include the candidates which are not in the Ind column of the origin pedigree when the cand parameter is not NULL.

# Changes in version 0.2.1 released on 17 Nov 2018
## New features
## Bug fixes
1. Repel the overlapping nodes due to very small differences (digits > 7) among x positions of nodes
