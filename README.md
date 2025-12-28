# visPedigree

Tidying and visualization for animal pedigrees.

The `visPedigree` package provides tools to check for duplicate and bisexual individuals, detect pedigree loops, add missing founders, sort parents before offspring, and trace the pedigrees of specified candidates. It generates hierarchical graphs for all individuals in a pedigree and can handle very large datasets (> 10,000 individuals per generation) by compacting full-sib groups. It is particularly effective for aquatic animal pedigrees, which often include numerous full-sib families per generation in nucleus breeding populations.

![simple pedigree](https://luansheng.netlify.app/post/2018-11-09-vispedigree-use-guide_files/figure-html/smallped-1.png)

More complex pedigree graphs can be found in the [vignettes](#vignettes).

## Installation

### From CRAN
*(Note: The package is preparing for submission and not yet available on CRAN)*
```R
# install.packages("visPedigree")
```

### From GitHub
You can install the development version from GitHub using the `devtools` package:
```R
# install.packages("devtools")
devtools::install_github("luansheng/visPedigree")
```

## Quick Start
```R   
library(visPedigree)

# Example 1: Tidy and visualize a small pedigree
# Draw the pedigree, compacting full-sib individuals and highlighting candidates.
# Note: the compacted candidates are also highlighted here.
cands <- c("Y", "Z1", "Z2")
small_ped |>
  tidyped(cand = cands) |>
  visped(compact = TRUE, highlight = cands)

# Example 2: Calculate and show inbreeding coefficients
library(data.table)
data.table(
  Ind = c("A", "B", "C", "D", "E"),
  Sire = c(NA, NA, "A", "C", "C"),
  Dam = c(NA, NA, "B", "B", "D"),
  Sex = c("male", "female", "male", "female", "male")
) |>
  tidyped(inbreed = TRUE) |>
  visped(highlight = c("E"), showf = TRUE)

```

## <a id="vignette">Vignette</a>
Drawing an animal pedigree using the visPedigree package [EN](https://luansheng.netlify.app/2018/11/09/vispedigree-use-guide/) [CN](https://luansheng.netlify.app/2018/09/24/the-first-package-vispedigree-0-1/)      

## Citation
LUAN Sheng (2018). visPedigree: A package for tidying and drawing animal pedigree. URL https://github.com/luansheng/visPedigree.

