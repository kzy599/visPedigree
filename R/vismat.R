#' Visualize Relationship Matrices
#'
#' @description \code{vismat} provides visualization tools for relationship 
#' matrices (A, D, AA, etc.). It supports individual-level heatmaps and 
#' relationship coefficient histograms.
#'
#' @param mat A relationship matrix (e.g., from \code{\link{pedmatrix}}), or a 
#' named list containing matrices, or a \code{tidyped} object.
#' @param ped Optional. A tidied pedigree object. Used for extracting labels 
#' or grouping information.
#' @param type Character, type of visualization. Supported:
#' \itemize{
#'   \item "heatmap": Individual or group-level relationship heatmap (default).
#'   \item "histogram": Histogram of relationship coefficients.
#' }
#' @param ids Character vector of individual IDs to display. If NULL, all 
#' individuals in the matrix are shown.
#' @param reorder Logical. If TRUE, rows and columns are reordered using 
#' hierarchical clustering to bring related individuals together. Only affects 
#' heatmap visualization. Skipped for large matrices (N > 2000).
#' @param grouping Optional. Name of the column in \code{ped} to group by 
#' (e.g., "Family", "Gen", "Year"). Used for grouping in heatmaps.
#' @param labelcex Numeric. Manual control for font size of individual 
#' labels.
#' @param ... Additional arguments passed to the plotting function 
#' (\code{lattice::levelplot} or \code{lattice::histogram}).
#'
#' @return No returned value. The plot is generated on the current device.
#' @export
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dist hclust
#' @importFrom lattice levelplot panel.levelplot panel.abline histogram
#' @importFrom data.table as.data.table
vismat <- function(mat, ped = NULL, type = "heatmap", ids = NULL, reorder = TRUE, grouping = NULL, labelcex = NULL, ...) {
  # 0. If input is a tidyped object, calculate A first
  if (inherits(mat, "tidyped")) {
    ped <- mat
    res <- pedmatrix(ped, method = "A")
    if (is.list(res) && "A" %in% names(res)) {
      mat <- res$A
    } else {
      mat <- res
    }
  }

  # 1. Handle list input from pedmatrix
  if (is.list(mat) && !is.matrix(mat) && !inherits(mat, "Matrix")) {
    preferred <- c("A", "D", "AA", "Ainv")
    found <- intersect(preferred, names(mat))
    if (length(found) > 0) {
      target_name <- found[1]
      if (length(mat) > 1) {
        message(sprintf("Found multiple matrices in list. Using '%s'.", target_name))
      }
      mat <- mat[[target_name]]
    } else {
      is_mat <- sapply(mat, function(x) is.matrix(x) || inherits(x, "Matrix"))
      if (any(is_mat)) {
        mat <- mat[[which(is_mat)[1]]]
      } else {
        found_types <- sapply(mat, function(x) class(x)[1])
        stop(sprintf("Input list does not contain any valid matrices. Found: %s", paste(found_types, collapse = ", ")))
      }
    }
  }

  # 2. Validation
  if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
    stop("Input 'mat' must be a matrix or a Matrix object.")
  }

  # 3. Subset by IDs if requested
  if (!is.null(ids)) {
    all_names <- rownames(mat)
    if (is.null(all_names)) {
      stop("Matrix has no rownames. Subset by 'ids' is not possible.")
    }
    valid_ids <- intersect(ids, all_names)
    if (length(valid_ids) == 0) {
      stop("None of the specified 'ids' were found in the matrix.")
    }
    mat <- mat[valid_ids, valid_ids, drop = FALSE]
  }

  # 3b. Reorder for clustering (if requested and type is heatmap)
  if (reorder && type == "heatmap") {
    n <- nrow(mat)
    if (n > 2000) {
      warning("Matrix too large for reordering (N > 2000). Skipping clustering.")
      reorder <- FALSE  # Mark as not reordered for downstream logic
    } else {
      # CRITICAL: We use dist() on the matrix rows (Relationship Profile Distance).
      # Siblings share almost identical relationship profiles with the whole population,
      # making their Euclidean distance much smaller than the Parent-Child distance.
      d_mat <- as.matrix(mat)
      dist_obj <- stats::dist(d_mat)
      hc <- stats::hclust(dist_obj, method = "ward.D2")
      mat <- mat[hc$order, hc$order, drop = FALSE]
    }
  }

  # 4. Handle Grouping (Aggregate)
  if (!is.null(grouping)) {
    if (is.null(ped)) {
      stop("'ped' must be provided when using 'grouping'.")
    }

    # Ensure ped is a data.table and includes the required column
    ped_dt <- data.table::as.data.table(ped)
    if (!grouping %in% names(ped_dt)) {
      stop(sprintf("Column '%s' not found in pedigree.", grouping))
    }

    # Match matrix IDs to groups
    # tidyped objects always use "Ind" as the ID column
    id_col <- "Ind"
    if (!id_col %in% names(ped_dt)) {
      # Fallback to the first column if "Ind" is missing (unlikely for tidyped)
      id_col <- names(ped_dt)[1]
    }
    mapping <- ped_dt[, .(id = as.character(get(id_col)), grp = get(grouping))]

    # Get group for each matrix row
    mat_ids <- rownames(mat)
    if (is.null(mat_ids)) stop("Matrix must have row names for grouping.")

    mat_grps <- mapping[match(mat_ids, id), grp]
    if (any(is.na(mat_grps))) {
      warning("Some individuals in matrix not found in pedigree. Using 'Unknown' group.")
      mat_grps[is.na(mat_grps)] <- "Unknown"
    }

    # Aggregate: Mean relationship between groups
    mat_grps <- as.character(mat_grps)
    grps_unique <- sort(unique(mat_grps))
    n_grp <- length(grps_unique)

    # Check if this is an inverse matrix (averaging inverse elements is usually not useful)
    # Note: Detection based on common naming patterns in row/col names or attributes
    mat_name <- if (!is.null(attr(mat, "method"))) attr(mat, "method") else ""
    is_inv <- grepl("inv", mat_name, ignore.case = TRUE) || 
              grepl("inv", paste(rownames(mat), collapse = ""), ignore.case = TRUE)
    if (is_inv) {
      warning("Aggregating an inverse relationship matrix. Mean values may not be biologically meaningful.")
    }

    message(sprintf("Aggregating %d individuals into %d groups based on '%s'...", nrow(mat), n_grp, grouping))

    agg_mat <- matrix(0, nrow = n_grp, ncol = n_grp, dimnames = list(grps_unique, grps_unique))

    # Vectorized aggregation using tapply-style indexing
    # Pre-compute all indices for efficiency
    idx_list <- lapply(grps_unique, function(g) which(mat_grps == g))
    
    # Efficient aggregation for both dense and sparse matrices
    for (i in seq_len(n_grp)) {
      idx_i <- idx_list[[i]]
      for (j in i:n_grp) {
        idx_j <- idx_list[[j]]
        # Extract sub-block and compute mean safely
        sub_block <- mat[idx_i, idx_j, drop = FALSE]
        val <- mean(as.matrix(sub_block), na.rm = TRUE)

        # Ensure value is finite
        if (!is.finite(val)) val <- 0

        agg_mat[i, j] <- val
        if (i != j) agg_mat[j, i] <- val  # Symmetry
      }
    }
    mat <- agg_mat
    # Note: User's reorder setting is preserved for grouped matrices

    # Set default titles for later use in heatmap
    grouping_main <- sprintf("Grouped Relationship Heatmap (%s)", grouping)
  } else {
    grouping_main <- NULL
  }

  # 5. Rendering
  if (type == "heatmap") {
    dots <- list(...)
    if (!"main" %in% names(dots)) {
      dots$main <- if(!is.null(grouping_main)) grouping_main else "Relationship Matrix Heatmap"
    }
    if (!"xlab" %in% names(dots)) dots$xlab <- if(!is.null(grouping)) grouping else (if(reorder) "Clustered Individuals" else "Individuals")
    if (!"ylab" %in% names(dots)) dots$ylab <- if(!is.null(grouping)) grouping else (if(reorder) "Clustered Individuals" else "Individuals")
    
    nature_genetics_palette <- grDevices::colorRampPalette(c(
      "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84", 
      "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"
    ))(100)
    if (!"col.regions" %in% names(dots)) {
      dots$col.regions <- nature_genetics_palette
    }

    # Add a color key (legend) by default
    if (!"colorkey" %in% names(dots)) {
      dots$colorkey <- TRUE
    }

    # IMPORTANT: Matrix::image() from the Matrix package handles dispatch in a way 
    # that often conflicts with lattice arguments for base matrices.
    # We use lattice::levelplot directly to ensure full grid (incl. 0s) and borders.
    mat_to_plot <- as.matrix(mat)
    
    # Final safety check for lattice
    if (!any(is.finite(mat_to_plot))) {
      stop("The relationship matrix contains no finite values to plot.")
    }
    
    # Set up scales for individual labels
    # We want Individual 1 at Top-Left, so Y-axis labels must be reversed.
    # We now allow up to 500 labels with dynamic font scaling.
    if (nrow(mat) <= 500 && !"scales" %in% names(dots)) {
      # Use manual labelcex if provided, otherwise dynamic font size
      use_cex <- if(!is.null(labelcex)) labelcex else max(0.2, 0.7 * (1 - (nrow(mat) - 20) / 600))
      
      dots$scales <- list(
        x = list(at = 1:nrow(mat), labels = rownames(mat), rot = 90, cex = use_cex),
        y = list(at = 1:ncol(mat), labels = rev(colnames(mat)), cex = use_cex)
      )
    }

    # Border/Grid logic: Use panel function to draw white grid lines.
    # Limit grid to N <= 100 to avoid visual clutter.
    if (nrow(mat) <= 100 && !"panel" %in% names(dots)) {
      dots$panel <- function(x, y, z, ...) {
        lattice::panel.levelplot(x, y, z, ...)
        lattice::panel.abline(h = (0:nrow(mat)) + 0.5, v = (0:nrow(mat)) + 0.5, col = "white", lwd = 0.4)
      }
    }

    # Execute plotting using lattice::levelplot
    # Transpose and reverse rows to make mat[1,1] appear at Top-Left (1, n)
    p <- do.call(lattice::levelplot, c(list(x = t(mat_to_plot[rev(seq_len(nrow(mat_to_plot))), , drop = FALSE])), dots))
    
    print(p)
    return(invisible(p))
    
  } else if (type == "histogram") {
    # Histogram of relationship coefficients
    # Useful for checking the range of inbreeding/kinship
    dots <- list(...)
    vals <- as.numeric(mat[lower.tri(mat)])
    
    if (!"main" %in% names(dots)) dots$main <- "Distribution of Relationship Coefficients (Kinship)"
    if (!"xlab" %in% names(dots)) dots$xlab <- "Relationship Value"
    if (!"ylab" %in% names(dots)) dots$ylab <- "Frequency (%)"
    
    p <- do.call(lattice::histogram, c(list(x = ~vals), dots))
    print(p)
    return(invisible(p))

  } else {
    stop(sprintf("Visualization type '%s' is not supported. Use 'heatmap' or 'histogram'.", type))
  }
}
