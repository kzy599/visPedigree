#' A large pedigree with big family sizes
#'
#' A dataset containing a pedigree with many full-sib individuals per family.
#'
#' @format A data.table with 7 columns:
#' \describe{
#'   \item{Ind}{Individual ID}
#'   \item{Sire}{Sire ID}
#'   \item{Dam}{Dam ID}
#'   \item{Sex}{Sex of the individual}
#'   \item{Year}{Year of birth}
#'   \item{IndNum}{Numeric ID for individual}
#'   \item{SireNum}{Numeric ID for sire}
#'   \item{DamNum}{Numeric ID for dam}
#' }
"big_family_size_ped"

#' A deep pedigree
#'
#' A dataset containing a pedigree with many generations.
#'
#' @format A data.table with 4 columns:
#' \describe{
#'   \item{Ind}{Individual ID}
#'   \item{Sire}{Sire ID}
#'   \item{Dam}{Dam ID}
#'   \item{Sex}{Sex of the individual}
#' }
"deep_ped"

#' A pedigree with loops
#'
#' A dataset containing a pedigree with circular mating loops.
#'
#' @format A data.table with 3 columns:
#' \describe{
#'   \item{Ind}{Individual ID}
#'   \item{Sire}{Sire ID}
#'   \item{Dam}{Dam ID}
#' }
"loop_ped"

#' A simple pedigree
#'
#' A small dataset containing a simple pedigree for demonstration.
#'
#' @format A data.table with 3 columns:
#' \describe{
#'   \item{Ind}{Individual ID}
#'   \item{Sire}{Sire ID}
#'   \item{Dam}{Dam ID}
#' }
"simple_ped"

#' A small pedigree
#'
#' A small dataset containing a pedigree with some missing parents.
#'
#' @format A data.frame with 3 columns:
#' \describe{
#'   \item{Ind}{Individual ID}
#'   \item{Sire}{Sire ID}
#'   \item{Dam}{Dam ID}
#' }
"small_ped"
