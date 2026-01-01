#' Internal constructor for tidyped class
#' @param x A data.table object
#' @return A tidyped object
#' @keywords internal
new_tidyped <- function(x) {
  stopifnot(inherits(x, "data.table"))
  attr(x, "tidyped") <- TRUE
  if (!inherits(x, "tidyped")) {
    class(x) <- c("tidyped", class(x))
  }
  x
}

#' Internal validator for tidyped class
#' @param x A tidyped object
#' @return The object if valid, otherwise an error
#' @keywords internal
validate_tidyped <- function(x) {
  if (!inherits(x, "tidyped")) {
    stop("Object must be of class 'tidyped'", call. = FALSE)
  }
  if (!inherits(x, "data.table")) {
    stop("Object must inherit from 'data.table'", call. = FALSE)
  }

  needed <- c("Ind", "Sire", "Dam", "Sex")
  missing_cols <- setdiff(needed, names(x))
  
  if (length(missing_cols) > 0) {
    stop(
      "The pedigree object is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ". Please ensure it has been processed by tidyped().",
      call. = FALSE
    )
  }
  
  if (!isTRUE(attr(x, "tidyped"))) {
    stop(
      "The object must have the 'tidyped' attribute set to TRUE.",
      call. = FALSE
    )
  }
  
  x
}

#' Print method for tidyped pedigree
#' @param x A tidyped object
#' @param ... Additional arguments (ignored)
#' @importFrom utils head
#' @export
print.tidyped <- function(x, ...) {
  cat("Tidy Pedigree Object\n")
  print(head(as.data.frame(x), 10))
  if (nrow(x) > 10) {
    cat("... with", nrow(x) - 10, "more rows\n")
  }
  invisible(x)
}

#' Plot a tidy pedigree
#'
#' @param x A \code{tidyped} object.
#' @param ... Additional arguments passed to \code{\link{visped}}.
#' @export
plot.tidyped <- function(x, ...) {
  validate_tidyped(x)
  visped(x, ...)
}

#' Summary method for tidyped objects
#'
#' @param object A tidyped object.
#' @param ... Additional arguments (ignored).
#' @return A summary.tidyped object containing pedigree statistics.
#' @export
summary.tidyped <- function(object, ...) {
  x <- object
  res <- list()
  
  res$n_ind <- nrow(x)
  res$n_male <- sum(x$Sex == "male", na.rm = TRUE)
  res$n_female <- sum(x$Sex == "female", na.rm = TRUE)
  res$n_unknown_sex <- sum(is.na(x$Sex))
  
  # Founders: Individuals with both parents unknown (NA)
  res$n_founders <- sum(is.na(x$Sire) & is.na(x$Dam))
  
  # Isolated: Gen == 0
  if ("Gen" %in% names(x)) {
    res$n_isolated <- sum(x$Gen == 0)
    res$max_gen <- max(x$Gen, na.rm = TRUE)
  } else {
    res$n_isolated <- NA
    res$max_gen <- NA
  }
  
  # Candidates
  if ("Cand" %in% names(x)) {
    res$n_cand <- sum(x$Cand, na.rm = TRUE)
  } else {
    res$n_cand <- NA
  }
  
  class(res) <- "summary.tidyped"
  res
}

#' Print method for summary.tidyped
#'
#' @param x A summary.tidyped object.
#' @param ... Additional arguments (ignored).
#' @export
print.summary.tidyped <- function(x, ...) {
  cat("Pedigree Summary:\n")
  cat("-----------------\n")
  cat("Total Individuals: ", x$n_ind, "\n")
  cat("  - Males:   ", x$n_male, "\n")
  cat("  - Females: ", x$n_female, "\n")
  if (x$n_unknown_sex > 0) {
    cat("  - Unknown Sex: ", x$n_unknown_sex, "\n")
  }
  cat("\n")
  cat("Founders (parents unknown): ", x$n_founders, "\n")
  
  if (!is.na(x$n_isolated) && x$n_isolated > 0) {
    cat("Isolated Individuals (Gen 0): ", x$n_isolated, "\n")
  }
  
  if (!is.na(x$max_gen)) {
    cat("Maximum Generation: ", x$max_gen, "\n")
  }
  
  if (!is.na(x$n_cand)) {
    cat("Candidates Traced: ", x$n_cand, "\n")
  }
  
  invisible(x)
}
