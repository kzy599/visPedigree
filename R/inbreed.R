#' Calculate inbreeding coefficients
#'
#' \code{inbreed} function calculates the inbreeding coefficients for all individuals in a tidied pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and calculates the inbreeding coefficients using optimized C++ code based on the Meuwissen & Luo (1992) algorithm. It is the core engine used by both \code{tidyped(..., inbreed = TRUE)} and \code{pedmatrix(..., method = "f")}, ensuring consistent results across the package. It is significantly faster than standard R implementations for large pedigrees.
#'
#' @param ped A \code{tidyped} object.
#' @param ... Additional arguments (currently ignored).
#' @return A \code{tidyped} object with an additional column \strong{f}.
#' @export
#' @import data.table
inbreed <- function(ped, ...) {
  validate_tidyped(ped)

  # Check for numeric pedigree columns
  if (!all(c("IndNum", "SireNum", "DamNum") %in% names(ped))) {
    ped <- tidyped(ped, addnum = TRUE)
  }

  # Ensure sorted by IndNum for the Trace algorithm
  ped_work <- copy(ped)
  setorder(ped_work, IndNum)

  # Call Rcpp function
  res_f <- cpp_calculate_inbreeding(ped_work$SireNum, ped_work$DamNum)
  
  # Map back to original ped
  ped_new <- copy(ped)
  # since we used ped_work order, map based on IndNum
  f_map <- data.table(IndNum = ped_work$IndNum, f_val = res_f$f)
  ped_new[f_map, f := i.f_val, on = "IndNum"]
  
  return(new_tidyped(ped_new))
}
