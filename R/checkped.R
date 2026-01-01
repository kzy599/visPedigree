#' Check a pedigree
#'
#' \code{checkped} function checks a pedigree.
#'
#' This function takes a pedigree, detects missing parents, checks for duplicated or bisexual individuals, adds missing founders, and sorts the pedigree. All individuals' sex will be inferred if there is no sexual information in the pedigree. If the pedigree includes the column \strong{Sex}, then individuals' sexes need to be recoded as "male", "female", or NA (unknown sex). Missing sexes will be identified from the pedigree structure and be added if possible.
#' @param ped A data.table or data frame with the pedigree, which includes the first three columns: \strong{individual}, \strong{sire} and \strong{dam} IDs. More columns, such as sex and generation can be included in the pedigree file. Names of the three columns can be assigned as you would like, but their orders must not be changed in the pedigree. Individual IDs should not be coded as "", " ", "0", "*", or "NA"; otherwise, these individuals will be removed from the pedigree. Missing parents should be denoted by either "NA", "0", or "*". Space and "" will also be recoded as missing parents, but are not recommended.
#' @param addgen A logical value indicating whether the generation number of an individual will be generated. The default value is TRUE, then a new column named \strong{Gen} will be added to the returned data.table.
#' @param ... Additional arguments passed to \code{\link{sortped}}.
#'
#' @return A data.table with the checked pedigree is returned. Individual, sire, and dam ID columns are renamed as \strong{Ind}, \strong{Sire} and \strong{Dam}. Missing parents are replaced with the default missing value \strong{NA}. The \strong{Sex} column contains the individual's sex (male, female, or NA for unknown). The \strong{Gen} column will be included when the parameter \emph{addgen} is TRUE. Ind, Sire, Dam, and Sex columns are character; The Gen column is an integer type.
#' @keywords internal
#' @import data.table
checkped <- function(ped, addgen = TRUE, ...) {
  ped_new <- copy(ped)
  setnames(ped_new,
           old = colnames(ped_new)[1:3],
           new = c("Ind", "Sire", "Dam"))
  setkey(ped_new, Ind, Sire, Dam)
  #===Detecting and setting missing values=============================================
  Ind <- Sire <- Dam <- NULL
  ped_new[, ":="(Ind = as.character(Ind),
                 Sire = as.character(Sire),
                 Dam = as.character(Dam))]
  # Individuals will be deleted if the Ind column contains "", " ", "0", "*", "NA", or NA.
  if (any(ped_new$Ind %in% c("", " ", "0", "*", "NA", NA))) {
    warning("Missing values in 'Ind' column; records discarded. Ensure 'Ind', 'Sire', 'Dam' are the first three columns.")
    ped_new <-
      ped_new[-which(ped_new$Ind %in% c("", " ", "0", "*", "NA", NA))]
  }
  # Missing parents are shown by "", " ", "0", "*", and "NA" in the pedigree file,
  # they are set as NA
  if (length(ped_new[Sire %in% c("", " ", "0", "*", "NA"), Sire]) > 0) {
    ped_new[Sire %in% c("", " ", "0", "*", "NA"), Sire := NA]
  }

  if (length(ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam]) > 0) {
    ped_new[Dam %in% c("", " ", "0", "*", "NA"), Dam := NA]
  }

  #The program will stop if no parents are in the sire and dam columns.
  if (all(is.na(ped_new$Sire)) & all(is.na(ped_new$Dam))) {
    stop("All parents are missing; no pedigree can be built. Please check your data.")
  }


  #====Detect or delete duplicated records===============================================
  # If the duplicated records have the same individual, sire, and dam ID,
  # one record will be kept, other records will be deleted.
  if (anyDuplicated(ped_new, by = c("Ind", "Sire", "Dam")) > 0) {
    ped_new_dup <-
      ped_new[duplicated(ped_new, by = c("Ind", "Sire", "Dam"))]
    ped_new_dup_num <- nrow(ped_new_dup)
    n = ped_new_dup_num
    if (n > 5) {
      n = 5
    }
    warning(sprintf("Removed %d records with duplicate Ind/Sire/Dam IDs. Showing first %d:", ped_new_dup_num, n))
    for (i in 1:n) {
      warning(paste(ped_new_dup[i], collapse = ", "))
    }
    if (ped_new_dup_num > 5) {
      warning("...")
    }

    ped_new <- unique(ped_new, by = c("Ind", "Sire", "Dam"))
  }
  # If the duplicated records only have the same individual ID,
  #  and their sire and dam IDs are different,
  # this program will stop because it is a fatal error.
  if (anyDuplicated(ped_new, by = c("Ind")) > 0) {
    ped_new_dup <-
      ped_new[duplicated(ped_new, by = c("Ind"))]
    ped_new_dup_num <- nrow(ped_new_dup)
    n = ped_new_dup_num
    if (n > 5) {
      n = 5
    }
    warning(sprintf("Found %d duplicate Ind IDs. Showing first %d:", ped_new_dup_num, n))

    for (i in 1:n) {
      warning(paste(ped_new_dup[i], collapse = ", "))
    }
    if (ped_new_dup_num > 5) {
      warning("...")
    }
    stop("Please check the pedigree!")
  }

  #===Find bisexual parents============================================================
  sires <- unique(ped_new$Sire)
  if (any(is.na(sires))) {
    sires <- sires[-which(is.na(sires))]
  }

  dams  <- unique(ped_new$Dam)
  if (any(is.na(dams))) {
    dams  <- dams[-which(is.na(dams))]
  }
  bisexual_parents <- sort(unique(sires[sires %chin% dams]))

  #===Renewing pedigree by adding missing parents or founders==========================
  sire_dam_vect <- unique(c(sires, dams))
  if (sum(!(sire_dam_vect %chin% ped_new$Ind)) > 0) {
    sire_dam_vect_missing <- sire_dam_vect[!(sire_dam_vect %chin% ped_new$Ind)]
    sire_dam_vect_missing_DT <- setDT(list(
      Ind = sire_dam_vect_missing,
      Sire = rep(NA, length(sire_dam_vect_missing)),
      Dam = rep(NA, length(sire_dam_vect_missing))
    ))
    ped_new <- rbind(sire_dam_vect_missing_DT, ped_new, fill = TRUE)
  }

  #===sorting parents in front of offspring in the individual column.
  SeqNumInd <- SeqNumSire <- SeqNumDam <- NULL
  ped_new[, SeqNumInd := .I]
  ped_new[,SeqNumSire:=SeqNumInd[match(Sire,Ind)]]
  ped_new[,SeqNumDam:=SeqNumInd[match(Dam,Ind)]]
  # Individuals are resorted if their orders are not right.
  if (any(
    c(
      ped_new$SeqNumInd < ped_new$SeqNumSire,
      ped_new$SeqNumInd < ped_new$SeqNumDam
    ),
    na.rm = TRUE
  ) | addgen)  {
    ped_new <- sortped(ped_new, addgen, ...)
  }

  attr(ped_new, "bisexual_parents") <- bisexual_parents

  # delete internal fields
  ped_new[,":="(SeqNumInd=NULL,SeqNumSire=NULL,SeqNumDam=NULL)]

  #===Add individual sex==========================================================
  col_names <- colnames(ped_new)
  Sex <- NULL
  if (!("Sex" %in% col_names)) {
    ped_new[, Sex := NA_character_]
    if (any(!is.na(ped_new$Sire))) {
      ped_new[Ind %chin% Sire, Sex:="male"]
    }
    if (any(!is.na(ped_new$Dam))) {
      ped_new[Ind %chin% Dam, Sex:="female"]
    }
  }
  if ("Sex" %in% col_names) {
    if (length(ped_new[Sex %in% c("", " ", "NA"), Sex]) > 0) {
      ped_new[Sex %in% c("", " ", "NA"), Sex := NA]
    }

    if (any(!is.na(ped_new$Sire))) {
      ped_new[is.na(Sex) & (Ind %chin% Sire), Sex:="male"]
    }
    if (any(!is.na(ped_new$Dam))) {
      ped_new[is.na(Sex) & (Ind %chin% Dam), Sex:="female"]
    }
  }
  if (any(c(!is.na(ped_new$Sire),!is.na(ped_new$Dam)))) {
    ped_new[!is.na(Sex),Sex:=tolower(Sex)]
    sex_name <- unique(ped_new[!is.na(Sex),Sex])
    if (!all(sex_name %in% c("male","female"))) {
      message("The Sex column contains values other than 'male' or 'female'!")
      message(paste("Sex:",paste(sex_name,collapse = " "),sep=" "))
    }
  }
  return(ped_new)
}
