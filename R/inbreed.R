#' Calculate inbreeding coefficients
#'
#' \code{inbreed} function calculates the inbreeding coefficients for all individuals in a tidied pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and calculates the inbreeding coefficients using the \code{makeDiiF} function from the \strong{nadiv} package. It prefers using numeric columns (\strong{IndNum}, \strong{SireNum}, \strong{DamNum}) if available, which is faster and more robust.
#'
#' @param ped A \code{tidyped} object.
#' @param ... Additional arguments passed to \code{\link[nadiv:makeDiiF]{makeDiiF}}.
#' @return A \code{tidyped} object with an additional column \strong{f}.
#' @export
#' @import data.table
#' @importFrom nadiv makeDiiF
inbreed <- function(ped, ...) {
  validate_tidyped(ped)

  if (!requireNamespace("nadiv", quietly = TRUE)) {
    stop("The 'nadiv' package is required for inbreeding calculations. Please install it.")
  }

  # Work on a copy and ensure it is sorted by IndNum (if available) 
  # to satisfy nadiv's requirement that parents appear before offspring.
  ped_work <- copy(ped)
  if ("IndNum" %in% colnames(ped_work)) {
    setorder(ped_work, IndNum)
  }

  # Prefer numeric IDs for nadiv if available (faster and more robust)
  if (all(c("IndNum", "SireNum", "DamNum") %in% colnames(ped_work))) {
    # Use numeric columns
    # nadiv expects ID, Dam, Sire
    nadiv_ped <- ped_work[, .(IndNum, DamNum, SireNum)]
    id_col <- "IndNum"
  } else {
    # Fallback to character IDs (nadiv will handle them)
    # nadiv expects ID, Dam, Sire
    nadiv_ped <- ped_work[, .(Ind, Dam, Sire)]
    id_col <- "Ind"
  }
  
  # Calculate inbreeding coefficients
  # Convert to data.frame to avoid data.table indexing issues in nadiv
  # nadiv::makeDiiF returns a list with 'f' vector matching the input ID column order
  # Suppress warnings from nadiv::numPed which redundantly warn about 0 as missing parent
  res <- suppressWarnings(nadiv::makeDiiF(as.data.frame(nadiv_ped), ...))
  
  # Create a mapping table to ensure correct order when joining back
  f_map <- data.table(id = nadiv_ped[[1]], f_val = res$f)
  
  ped_new <- copy(ped)
  # Join by ID to ensure f values are assigned to the correct individuals
  # regardless of the row order in the input 'ped'
  if (id_col == "IndNum") {
    ped_new[f_map, f := i.f_val, on = .(IndNum = id)]
  } else {
    ped_new[f_map, f := i.f_val, on = .(Ind = id)]
  }
  
  return(new_tidyped(ped_new))
}
