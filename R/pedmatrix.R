#' Genetic Relationship Matrices and Inbreeding Coefficients
#'
#' @description Optimized calculation of additive (A), its inverse (Ainv), 
#' dominance (D) matrices, and inbreeding coefficients (f) using Rcpp and 
#' Meuwissen & Luo (1992) algorithm.
#'
#' @param ped A tidied pedigree (data.table or data.frame).
#' @param method Character, one of:
#' \itemize{
#'   \item "f": Inbreeding coefficients
#'   \item "A": Additive relationship matrix
#'   \item "Ainv": Inverse of A (Henderson's rules)
#'   \item "D": Dominance relationship matrix
#'   \item "Dinv": Inverse of D (see Performance Notes below)
#'   \item "AA": Additive-by-additive relationship matrix (A # A)
#'   \item "AAinv": Inverse of AA
#' }
#' 
#' @details 
#' 
#' \strong{General Usage:}
#' 
#' Only a single method may be requested per call. This design prevents
#' accidental heavy computations and makes the API more explicit. If multiple
#' results are needed, call \code{pedmatrix()} separately for each method.
#'
#' The function returns the requested matrix or vector directly (not wrapped
#' in a list).
#' 
#' \strong{Performance Notes:}
#' 
#' \itemize{
#'   \item \strong{Ainv}: O(n) complexity using Henderson's direct rules. 
#'         Highly efficient for all pedigree sizes.
#'   \item \strong{Dinv/AAinv}: O(n³) complexity using dense matrix inversion. 
#'         Performance characteristics:
#'         \itemize{
#'           \item n < 500: ~10-20 ms (fast, recommended)
#'           \item n = 1000: ~40-60 ms (acceptable)
#'           \item n = 2000: ~130-150 ms (with Cholesky optimization)
#'           \item n > 2000: May become slow and memory intensive
#'         }
#'   \item \strong{Optimization}: By default (invert_method="auto"), the function
#'         automatically detects symmetric positive-definite matrices and uses
#'         Cholesky decomposition, which is ~2x faster than general LU decomposition.
#'         D and AA matrices are always symmetric positive-definite, so this
#'         optimization is always applied.
#'   \item \strong{Note}: Unlike Ainv which uses direct construction rules,
#'         Dinv/AAinv requires computing the matrix first and then inverting it. 
#'         For large pedigrees (n > 2000), consider whether the inverse is truly 
#'         needed or if alternative formulations are possible.
#' }
#' 
#' @param sparse Logical, if TRUE returns sparse Matrix (for Ainv, A, D, AA).
#' @param invert_method Character, method for matrix inversion (Dinv/AAinv):
#' \itemize{
#'   \item "auto": Auto-detect matrix properties and use optimal method (default)
#'   \item "sympd": Force Cholesky decomposition (symmetric positive-definite)
#'   \item "general": Force general LU decomposition
#' }
#' @param n_threads Integer, number of threads for parallel tasks. Currently
#' unused; C++ implementation uses all available cores automatically.
#'
#' @return Returns the requested object directly:
#' \itemize{
#'   \item For "f": Named numeric vector of inbreeding coefficients
#'   \item For "A", "D", "AA": Matrix (sparse or dense based on \code{sparse})
#'   \item For "Ainv", "Dinv", "AAinv": Matrix (sparse or dense)
#' }
#' @export
#'
#' @importFrom Matrix sparseMatrix diag
#' @useDynLib visPedigree, .registration = TRUE
#' @importFrom Rcpp evalCpp
pedmatrix <- function(ped, method = "f", sparse = TRUE, invert_method = "auto", n_threads = 0) {
  # Coerce method to character and validate single method behavior
  if (length(method) != 1) {
    stop("Only a single 'method' may be requested. Please call pedmatrix() separately for each desired method.")
  }
  method <- as.character(method)
  
  # Validate invert_method
  invert_method <- match.arg(invert_method, c("auto", "sympd", "general"))

  # 1. Standardize pedigree
  if (!inherits(ped, "tidyped")) {
    ped <- tidyped(ped, addnum = TRUE)
  }
  
  # Ensure sorted by IndNum to match matrix indexing
  ped <- data.table::copy(ped)
  data.table::setorder(ped, IndNum)

  sire <- ped$SireNum
  dam <- ped$DamNum
  n <- nrow(ped)
  
  # Use a named list to store results
  output <- list()
  
  # Lazy-loading objects to avoid redundant calculations
  # Note: A_dense, D_dense, etc. will fail for 1M individuals. 
  # We only build them if specifically requested in 'method'.
  cache <- list(f_res = NULL, A_dense = NULL, D_dense = NULL, AA_dense = NULL)
  
  get_f_res <- function() {
    if (is.null(cache$f_res)) {
      # Reverting to highly optimized sequential version for f
      # It is very fast and has no parallel overhead or sorting overhead
      cache$f_res <<- cpp_calculate_inbreeding(sire, dam)
    }
    cache$f_res
  }
  
  get_A_dense <- function() {
    if (n > 25000) stop("Pedigree too large (>25k) for dense A matrix. Use 'Ainv' or sparse methods.")
    if (is.null(cache$A_dense)) cache$A_dense <<- cpp_calculate_A(sire, dam)
    cache$A_dense
  }

  get_D_dense <- function() {
    if (n > 25000) stop("Pedigree too large (>25k) for dense D matrix.")
    if (is.null(cache$D_dense)) {
      A_mat <- get_A_dense()
      cache$D_dense <<- cpp_calculate_D(sire, dam, A_mat)
    }
    cache$D_dense
  }

  get_AA_dense <- function() {
    if (n > 25000) stop("Pedigree too large (>25k) for dense AA matrix.")
    if (is.null(cache$AA_dense)) {
      A_mat <- get_A_dense()
      cache$AA_dense <<- cpp_calculate_AA(A_mat)
    }
    cache$AA_dense
  }

  # Process methods (ordered by dependency)
  all_methods <- c("f", "A", "Ainv", "D", "Dinv", "AA", "AAinv")
  active_methods <- intersect(all_methods, method)

  # Check for invalid methods
  invalid <- setdiff(method, all_methods)
  if (length(invalid) > 0) {
    stop(sprintf("Invalid method(s): %s", paste(invalid, collapse = ", ")))
  }

  # Process the single requested method (previously supported multiple)
  for (m in active_methods) {
    if (m == "f") {
      # Prefer getting f from existing A diag, else calculated
      if (!is.null(cache$A_dense)) {
        output$f <- diag(cache$A_dense) - 1.0
      } else {
        output$f <- get_f_res()$f
      }
      
    } else if (m == "A") {
      A_mat <- get_A_dense()
      output$A <- if (sparse) Matrix::Matrix(A_mat, sparse = TRUE) else A_mat
      
    } else if (m == "D") {
      D_mat <- get_D_dense()
      output$D <- if (sparse) Matrix::Matrix(D_mat, sparse = TRUE) else D_mat
      
    } else if (m == "Dinv") {
      D_mat <- get_D_dense()
      Dinv <- switch(invert_method,
        "auto" = cpp_invert_auto(D_mat),
        "sympd" = cpp_invert_sympd(D_mat),
        "general" = cpp_invert_dense(D_mat)
      )
      output$Dinv <- if (sparse) Matrix::Matrix(Dinv, sparse = TRUE) else Dinv

    } else if (m == "AA") {
      AA_mat <- get_AA_dense()
      output$AA <- if (sparse) Matrix::Matrix(AA_mat, sparse = TRUE) else AA_mat

    } else if (m == "AAinv") {
      AA_mat <- get_AA_dense()
      AAinv <- switch(invert_method,
        "auto" = cpp_invert_auto(AA_mat),
        "sympd" = cpp_invert_sympd(AA_mat),
        "general" = cpp_invert_dense(AA_mat)
      )
      output$AAinv <- if (sparse) Matrix::Matrix(AAinv, sparse = TRUE) else AAinv
      
    } else if (m == "Ainv") {
      f_res <- get_f_res()
      # Pass n_threads if we want, but currently C++ uses max
      triplets <- cpp_build_ainv_triplets(sire, dam, f_res$dii)
      Ainv <- Matrix::sparseMatrix(
        i = triplets$i, j = triplets$j, x = triplets$v,
        dims = c(n, n), symmetric = TRUE
      )
      output$Ainv <- if (sparse) Ainv else as.matrix(Ainv)
    }
  }

  # 3. Add dimension names to matrices
  ind_names <- ped$Ind
  for (m in names(output)) {
    obj <- output[[m]]
    if (is.matrix(obj) || inherits(obj, "Matrix")) {
      # Square matrices: set both dimensions
      rownames(obj) <- ind_names
      colnames(obj) <- ind_names
      output[[m]] <- obj
    } else if (is.numeric(obj) && length(obj) == n) {
      # Vectors (f)
      names(obj) <- ind_names
      output[[m]] <- obj
    }
  }

  # If a single method was requested, simplify return to direct object
  if (length(active_methods) == 1) {
    return(output[[active_methods]])
  }

  return(output)
}
