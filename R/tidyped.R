#' Tidy and prepare a pedigree
#'
#' The \code{tidyped} function checks, sorts, traces, and returns a trimmed pedigree.
#'
#' This function takes a pedigree, checks for duplicate and bisexual individuals, detects pedigree loops, adds missing founders, sorts the pedigree, and traces the pedigree of specified candidates. If the \code{cand} parameter contains individual IDs, only those individuals and their ancestors or descendants are retained. Tracing direction and the number of generations can be specified using the \code{trace} and \code{tracegen} parameters. Virtual generations are inferred and assigned when \code{addgen = TRUE}. A numeric pedigree is generated when \code{addnum = TRUE}. Sex will be inferred for all individuals if sex information is missing. If the pedigree includes a \strong{Sex} column, values should be coded as "male", "female", or NA (unknown). Missing sex information will be inferred from the pedigree structure where possible.
#'
#' @param ped A data.table or data frame containing the pedigree. The first three columns must be \strong{individual}, \strong{sire}, and \strong{dam} IDs. Additional columns, such as sex or generation, can be included. Column names can be customized, but their order must remain unchanged. Individual IDs should not be coded as "", " ", "0", "*", or "NA"; otherwise, they will be removed. Missing parents should be denoted by "NA", "0", or "*". Spaces and empty strings ("") are also treated as missing parents but are not recommended.
#' @param cand A character vector of individual IDs, or NULL. If provided, only the candidates and their ancestors/descendants are retained.
#' @param trace A character value specifying the tracing direction: "\strong{up}", "\strong{down}", or "\strong{all}". "up" traces ancestors; "down" traces descendants; "all" traces both simultaneously. This parameter is only used if \code{cand} is not NULL. Default is "up".
#' @param tracegen An integer specifying the number of generations to trace. This parameter is only used if \code{trace} is not NULL. If NULL or 0, all available generations are traced.
#' @param addgen A logical value indicating whether to generate generation numbers. Default is TRUE, which adds a \strong{Gen} column to the output.
#' @param addnum A logical value indicating whether to generate a numeric pedigree. Default is TRUE, which adds \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns to the output.
#' @param inbreed A logical value indicating whether to calculate inbreeding coefficients. Default is FALSE. If TRUE, an \strong{f} column is added to the output.
#' @param ... Additional arguments passed to internal methods such as \code{\link{checkped}}, \code{\link{numped}}, and \code{\link{inbreed}}.
#'
#' @return A \code{tidyped} object (which inherits from \code{data.table}). Individual, sire, and dam ID columns are renamed to \strong{Ind}, \strong{Sire}, and \strong{Dam}. Missing parents are replaced with \strong{NA}. The \strong{Sex} column contains "male", "female", or NA. The \strong{Cand} column is included if \code{cand} is not NULL. The \strong{Gen} column is included if \code{addgen} is TRUE. The \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns are included if \code{addnum} is TRUE. The \strong{f} column is included if \code{inbreed} is TRUE. \strong{Ind}, \strong{Sire}, \strong{Dam}, and \strong{Sex} are character columns; \strong{Cand} is logical; \strong{Gen}, \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} are integer.
#' @seealso \code{\link{summary.tidyped}} for summarizing the pedigree.
#'
#' @examples
#' library(visPedigree)
#' library(data.table)
#'
#' # Tidy a simple pedigree
#' tidy_ped <- tidyped(simple_ped)
#' head(tidy_ped)
#'
#' # Calculate inbreeding coefficients
#' tidy_ped_inbreed <- tidyped(simple_ped, inbreed = TRUE)
#' head(tidy_ped_inbreed)
#'
#' # Trace ancestors of a specific individual
#' tidy_ped_J5X804 <- tidyped(simple_ped, cand = "J5X804")
#' head(tidy_ped_J5X804)
#'
#' # Trace ancestors back two generations
#' tidy_ped_J5X804_up_2 <- tidyped(simple_ped, cand = "J5X804", trace = "up", tracegen = 2)
#' tidy_ped_J5X804_up_2
#'
#' # Trace descendants of a specific individual
#' tidy_ped_J3Y620_down <- tidyped(simple_ped, cand = "J3Y620", trace = "down")
#' tidy_ped_J3Y620_down
#'
#' # Trace descendants down two generations
#' tidy_ped_J3Y620_down_2 <- tidyped(simple_ped, cand = "J3Y620", trace = "down", tracegen = 2)
#' tidy_ped_J3Y620_down_2
#'
#' @import data.table
#' @export
tidyped <-
  function(ped,
           cand = NULL,
           trace = "up",
           tracegen = NULL,
           addgen = TRUE,
           addnum = TRUE,
           inbreed = FALSE,
           ...) {
    ped_is_DT <- "data.table" %in% class(ped)
    if (!ped_is_DT) {
      ped_inter <- as.data.table(ped)
    } else {
      ped_inter <- copy(ped)
    }
    attr(ped_inter, "tidyped") <- FALSE

    if (!is.null(cand)) {
      if (is.character(trace)) {
        if (!trace %in% c("up", "down", "all")) {
          stop("The trace parameter only is assigned using \"up\", \"down\" or \"all\"!")
        }
      } else {
        stop("The trace parameter only is assigned using \"up\", \"down\" or \"all\"!")
      }
    }

    if (!is.null(tracegen)) {
      if (is.numeric(tracegen)) {
        if (!tracegen == floor(tracegen)) {
          stop("The tracegen parameter need to be an integer!")
        }
      } else {
        stop("The tracegen parameter need to be an integer!")
      }
    }

    if (!is.logical(addgen)) {
      stop("The addgen parameter only is assigned using TRUE or FALSE")
    }

    if (!is.logical(addnum)) {
      stop("The addnum parameter only is assigned using TRUE or FALSE")
    }

    ped_colnames <- colnames(ped_inter)

    # Delete Cand column
    if (!is.null(cand)) {
      if (c("Cand") %in% ped_colnames) {
        ped_inter[, Cand := NULL]
        warning("The column Cand of the original pedigree is deleted.")
      }
    }

    # The Gen column will be deleted
    if (c("Gen") %in% ped_colnames) {
      ped_inter[, Gen := NULL]
      warning("The column Gen of the original pedigree is deleted.")
    }

    # IndNum SireNum, DamNum or f columns will be deleted if they are to be re-generated
    cols_to_delete <- c()
    if (addnum) {
      cols_to_delete <- c(cols_to_delete, "IndNum", "SireNum", "DamNum")
    }
    if (inbreed) {
      cols_to_delete <- c(cols_to_delete, "f")
    }

    if (length(cols_to_delete) > 0 && any(cols_to_delete %in% ped_colnames)) {
      exist_columns <- cols_to_delete[cols_to_delete %in% ped_colnames]
      ped_inter[, (exist_columns) := NULL]
      warning(
        "The columns ",
        paste(exist_columns, collapse = ","),
        " of the original pedigree are deleted."
      )
    }

    # Check pedigree and sort it
    ped_inter <- withCallingHandlers(
      checkped(ped_inter, addgen = addgen, ...),
      warning = function(w) {
        if (grepl("shallow copy of this data.table", conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
    bisexual_parents <- attr(ped_inter, "bisexual_parents")
    if (!is.null(bisexual_parents) && length(bisexual_parents) > 0) {
      warning(
        paste(
          "bisexual individuals (found in both Sire and Dam columns):",
          paste(bisexual_parents, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    # Pruning the pedigree by candidate
    if (!is.null(cand)) {
      if (!all(class(cand) %in% c("character"))) {
        cand <- as.character(cand)
      }

      if (!any(cand %in% ped_inter$Ind)) {
        stop("No candidates are in the pedigree!")
      } else {
        ped_inter[Ind %chin% cand, Cand := TRUE]
        ped_inter[!(Ind %chin% cand), Cand := FALSE]
      }
      
      # Use numped internally for tracing
      ped_num <- numped(ped_inter, ...)
      cand_num <- match(cand, ped_num$Ind, nomatch = 0)
      
      keep_ind <- c()
      
      if (trace %in% c("up", "all")) {
        # Trace from candidate to ancestry
        i <- 1
        keep_ind_backward <- cand_num[cand_num > 0]
        ind_n <- length(keep_ind_backward) + 1
        while (length(keep_ind_backward) != ind_n) {
          if (!is.null(tracegen) && i > tracegen) break
          ind_n <- length(keep_ind_backward)
          keep_ind_backward <- unique(c(unlist(ped_num[keep_ind_backward, .(SireNum, DamNum)]), keep_ind_backward))
          keep_ind_backward <- keep_ind_backward[keep_ind_backward > 0]
          i <- i + 1
        }
        keep_ind <- unique(c(keep_ind, keep_ind_backward))
      }

      if (trace %in% c("down", "all")) {
        # Trace from candidate to descendant
        i <- 1
        keep_ind_foreward <- cand_num[cand_num > 0]
        ind_n <- length(keep_ind_foreward) + 1
        while (length(keep_ind_foreward) != ind_n) {
          if (!is.null(tracegen) && i > tracegen) break
          ind_n <- length(keep_ind_foreward)
          keep_ind_foreward <- unique(c(ped_num[SireNum %in% keep_ind_foreward | DamNum %in% keep_ind_foreward, IndNum], keep_ind_foreward))
          keep_ind_foreward <- keep_ind_foreward[keep_ind_foreward > 0]
          i <- i + 1
        }
        keep_ind <- unique(c(keep_ind, keep_ind_foreward))
      }

      ped_new <- ped_num[sort(keep_ind)]
      
      # Ensure parents not in the traced set are set to NA to respect tracegen depth
      ped_new[!(SireNum %in% keep_ind), Sire := NA]
      ped_new[!(DamNum %in% keep_ind), Dam := NA]

      ped_new[, ":="(IndNum = NULL, SireNum = NULL, DamNum = NULL)]
      
      # Re-check to ensure consistency after pruning
      ped_inter <- checkped(ped_new, addgen = addgen, ...)
    }

    # Converting individual, sire, and dam IDs to sequential integer
    if (addnum) {
      ped_inter <- numped(ped_inter, ...)
    }

    # Add a new column Cand if it doesn't exist
    if (!is.null(cand) && !"Cand" %in% colnames(ped_inter)) {
      ped_inter[, Cand := Ind %in% cand]
    }
    
    # Convert to tidyped object
    ped_inter <- new_tidyped(ped_inter)

    # Calculate inbreeding coefficients using inbreed() function
    if (inbreed) {
      ped_inter <- inbreed(ped_inter, ...)
    }

    return(ped_inter)
  }

