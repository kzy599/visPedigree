# Changes in version 0.4 released on 24 Dec 2025
## New features
1. Added `inbreed` parameter to `tidyped()` function to calculate inbreeding coefficients using the `nadiv` package. The result is stored in a concise column named `f`.
2. Added `showf` parameter to `visped()` function to display inbreeding coefficients on the pedigree graph. This parameter only controls display; it does not trigger new calculations, ensuring better performance and separation of concerns.
3. Refactored `inbreed()` function as a standalone tool that operates on `tidyped` objects, utilizing numeric IDs for faster and more robust calculation via `nadiv::makeDiiF`.
4. Optimized `repeloverlap()` function using `data.table` for significantly better performance and 100% functional consistency with previous versions.
5. Added `highlight` parameter to `visped()` function. Users can now highlight specific individuals using a character vector of IDs or a list for custom colors.

## Bug fixes
1. Fixed spelling errors in `visped()` documentation ("shwon" -> "shown", "genertion" -> "generation").

// ...existing code...

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
