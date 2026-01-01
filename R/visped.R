#' Visualize a tidy pedigree
#'
#' \code{visped} function draws a graph of a full or compact pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function and outputs a hierarchical graph for all individuals in the pedigree. The graph can be shown on the default graphic device or saved as a PDF file. The PDF output is a vector drawing that is legible and avoids overlapping labels. It is especially useful when the number of individuals is large and individual labels are long. This function can draw the graph of a very large pedigree (> 10,000 individuals per generation) by compacting full-sib individuals. It is highly effective for aquatic animal pedigrees, which usually include many full-sib families per generation in nucleus breeding populations. The outline of a pedigree without individual labels is still shown if the width of a pedigree graph exceeds the maximum width (200 inches) of the PDF file.
#'
#' In the graph, two shapes and three colors are used. Circles represent individuals, and squares represent families. Dark sky blue indicates males, dark goldenrod indicates females, and dark olive green indicates unknown sex. For example, a dark sky blue circle represents a male individual; a dark goldenrod square represents all female individuals in a full-sib family when \code{compact = TRUE}.
#'
#' @param ped A \code{tidyped} object (which inherits from \code{data.table}). It is recommended that the pedigree is tidied and pruned by candidates using the \code{\link{tidyped}} function with the non-null parameter \code{cand}.
#' @param compact A logical value indicating whether IDs of full-sib individuals in one generation will be removed and replaced with the number of full-sib individuals. For example, if there are 100 full-sib individuals in one generation, they will be replaced with a single label "100" when \code{compact = TRUE}. The default value is FALSE.
#' @param outline A logical value indicating whether shapes without labels will be shown. A graph of the pedigree without individual labels is shown when setting \code{outline = TRUE}. This is useful for viewing the pedigree outline and identifying immigrant individuals in each generation when the graph width exceeds the maximum PDF width (200 inches). The default value is FALSE.
#' @param cex NULL or a numeric value changing the size of individual labels shown in the graph. \emph{cex} is an abbreviation for 'character expansion factor'. The \code{visped} function will attempt to estimate (\code{cex=NULL}) the appropriate cex value and report it in the messages. Based on the reported cex from a previous run, this parameter should be increased if labels are wider than their shapes in the PDF; conversely, it should be decreased if labels are narrower than their shapes. The default value is NULL.
#' @param showgraph A logical value indicating whether a plot will be shown in the default graphic device (e.g., the Plots panel in RStudio). This is useful for quick viewing without opening a PDF file. However, the graph on the default device may not be legible (e.g., overlapping labels or aliasing lines) due to size restrictions. It is recommended to set \code{showgraph = FALSE} for large pedigrees. The default value is TRUE.
#' @param file NULL or a character value specifying whether the pedigree graph will be saved as a PDF file. The PDF output is a legible vector drawing where labels do not overlap, even with many individuals or long labels. It is recommended to save the pedigree graph as a PDF file. The default value is NULL.
#' @param highlight NULL, a character vector of individual IDs, or a list specifying individuals to highlight. If a character vector is provided, individuals will be highlighted with a purple border while preserving their sex-based fill color. If a list is provided, it should contain:
#' \itemize{
#'   \item \code{ids}: (required) character vector of individual IDs to highlight.
#'   \item \code{frame.color}: (optional) hex color for the border of focal individuals.
#'   \item \code{color}: (optional) hex color for the fill of focal individuals.
#'   \item \code{rel.frame.color}: (optional) hex color for the border of relatives (used when \code{trace} is not NULL).
#'   \item \code{rel.color}: (optional) hex color for the fill of relatives (used when \code{trace} is not NULL).
#' }
#' For example: \code{c("A", "B")} or \code{list(ids = c("A", "B"), frame.color = "#9c27b0")}. The function will check if the specified individuals exist in the pedigree and issue a warning for any missing IDs. The default value is NULL.
#' @param trace A logical value or a character string. If TRUE, all ancestors and descendants of the individuals specified in \code{highlight} will be highlighted. If a character string, it specifies the tracing direction: "\strong{up}" (ancestors), "\strong{down}" (descendants), or "\strong{all}" (both). This is useful for focusing on specific families within a large pedigree. The default value is FALSE.
#' @param showf A logical value indicating whether inbreeding coefficients will be shown in the graph. If \code{showf = TRUE} and the column \strong{f} exists in the pedigree, the inbreeding coefficient will be appended to the individual label, e.g., "ID (0.05)". The default value is FALSE.
#' @param ... Additional arguments passed to \code{\link[igraph:plot.igraph]{plot.igraph}}.
#' @return No returned values. The graph will be plotted directly on graphic devices.
#' @note Isolated individuals (those with no parents and no progeny, assigned Gen 0) are automatically filtered out and not shown in the plot. A warning will be issued if any such individuals are removed.
#'
#' @examples
#' library(data.table)
#' # Drawing a simple pedigree
#' simple_ped
#' simple_ped_tidy <- tidyped(simple_ped)
#' visped(simple_ped_tidy)
#' # Highlighting an individual and its relatives
#' visped(simple_ped_tidy, highlight = "J5X804", trace = TRUE)
#' # Showing inbreeding coefficients in the graph
#' visped(simple_ped_tidy, showf = TRUE)
#' # Drawing a simple pedigree of a individual with id of J5X804
#' simple_ped_J5X804_tidy <- tidyped(simple_ped,cand=c("J5X804"))
#' visped(simple_ped_J5X804_tidy)
#' # Drawing the graph in the pdf file
#' visped(simple_ped_J5X804_tidy, file = tempfile(fileext = ".pdf"))
#' # Highlighting specific individuals with default colors
#' visped(simple_ped_tidy, highlight = c("Y", "Z1", "Z2"))
#' # Highlighting specific individuals with custom colors
#' visped(simple_ped_tidy, 
#'        highlight = list(ids = c("Y", "Z1"), frame.color = "#4caf50", color = "#81c784"))
#' # Drawing a compact pedigree
#' # The candidates' labels in 2007
#' cand_labels <- big_family_size_ped[(Year == 2007) & (substr(Ind,1,2) == "G8"),Ind]
#' big_ped_tidy <- tidyped(big_family_size_ped, cand = cand_labels)
#' \donttest{
#' visped(big_ped_tidy, compact = TRUE)
#' visped(big_ped_tidy, compact = TRUE, file = tempfile(fileext = ".pdf"))
#' # Individual labels are not shown
#' visped(big_ped_tidy, compact = TRUE, outline = TRUE, file = tempfile(fileext = ".pdf"))
#' }
#'
#' @import data.table
#' @import igraph
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics strwidth
#' @export
visped <- function(ped,
                   compact = FALSE, outline = FALSE, cex = NULL, showgraph = TRUE, file = NULL, 
                   highlight = NULL, trace = FALSE, showf = FALSE, ...) {
  validate_tidyped(ped)

  if (showf && !"f" %in% colnames(ped)) {
    warning("Inbreeding coefficients ('f' column) not found in pedigree. Please run tidyped(..., inbreed = TRUE) to calculate them.")
    showf <- FALSE
  }

  if (!isFALSE(trace)) {
    if (isTRUE(trace)) {
      trace <- "all"
    } else if (!trace %in% c("up", "down", "all")) {
      stop("trace must be one of FALSE, TRUE, 'up', 'down', or 'all'.")
    }
  }

  # Prepare graph data
  graph_data <- prepare_ped_graph(
    ped = ped,
    compact = compact,
    outline = outline,
    cex = cex,
    highlight = highlight,
    trace = trace,
    showf = showf,
    ...
  )

  g <- graph_data$g
  l <- graph_data$layout
  node_size <- graph_data$node_size

  #===Draw the pedigree================================================================
  if (showgraph) {
    plot_ped_igraph(g, l, node_size, ...)
  }

  if (!is.null(file)) {
    pdf(file = file,
        width = graph_data$canvas_width,
        height = graph_data$canvas_height)
    # Ensure the device is closed even if plotting fails.
    on.exit(dev.off(), add = TRUE)
    plot_ped_igraph(g, l, node_size, ...)
    message(paste("Pedigree saved to: ", file.path(getwd(), file), sep = ""))
  }

  if ((showgraph || !is.null(file)) && !outline) {
    current_cex <- if (is.null(cex)) graph_data$best_cex else cex
    message(paste("Label cex: ", current_cex, ". Adjust 'cex' if labels are too large or small.", sep = ""))
  }

  if ((showgraph || !is.null(file)) && is.null(file)) {
    message("Tip: Use 'file' to save as a legible vector PDF.")
  }

  if (showf && "f" %in% colnames(ped) && any(ped$f == 0, na.rm = TRUE)) {
    message("Note: Inbreeding coefficients of 0 are not shown in the graph.")
  }

  invisible(graph_data)
}

#' Prepare pedigree graph data
#'
#' Internal function to calculate coordinates and layout for pedigree visualization.
#'
#' @inheritParams visped
#' @return A list containing the igraph object, layout, and canvas dimensions.
#' @keywords internal
prepare_ped_graph <- function(ped, compact = FALSE, outline = FALSE, cex = NULL, 
                              highlight = NULL, trace = FALSE, showf = FALSE, ...) {
  ped_new <- copy(ped)

  # Ensure required columns for plotting exist
  if (!"Gen" %in% colnames(ped_new)) {
    ped_new <- sortped(ped_new, addgen = TRUE, ...)
  }
  if (!all(c("IndNum", "SireNum", "DamNum") %in% colnames(ped_new))) {
    ped_new <- numped(ped_new, ...)
  }

  # Reserved digits
  fixed_digits <- 7
  # Digits when calculating
  old_digits <- getOption("digits")
  options(digits = 20)
  on.exit(options(digits = old_digits), add = TRUE)

  ped_igraph_data <- ped2igraph(ped_new, compact, highlight, trace, showf)
  ped_igraph <- ped_igraph_data
  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  
  # Filter out isolated individuals (Gen == 0)
  if (any(real_node$gen == 0)) {
    n_isolated <- sum(real_node$gen == 0)
    warning(sprintf("Removed %d isolated individuals (no parents, no progeny) from the plot.", n_isolated))
    
    # Remove from real_node
    real_node <- real_node[gen != 0]
    
    # Remove from ped_igraph$node
    ped_igraph$node <- ped_igraph$node[gen != 0]
    
    # Remove from ped_igraph$edge
    valid_ids <- ped_igraph$node$id
    ped_igraph$edge <- ped_igraph$edge[from %in% valid_ids & to %in% valid_ids]
  }
  
  if (nrow(real_node) == 0) {
    stop("Pedigree contains no individuals to plot.")
  }

  gen_node_num <- real_node[, .N, by = gen]
  gen_max_size <- max(gen_node_num$N, na.rm = TRUE)

  #=== Obtaining the maximum width of nodes' label =====================================================
  # The maximum width of a PDF file is 200 inch
  pdf_max_width = 200  # Max PDF width in inches
  cexs <-  seq(from = 0.1, to = 1, by = 0.05)  # CEX sequence
  best_cex <- 0
  max_strwidth_label <- real_node[which.max(strwidth(real_node$label, cex = 1, units = "inches")), label]
  
  for (i in length(cexs):1) {
    # Obtaining the maximum width of a node's label: inch
    label_max_width <- max(strwidth(max_strwidth_label, cex = cexs[i], units = "inches"),
          na.rm = TRUE)
    # Fixing the width of the node when the number of nodes (individuals) is small in one generation
    # The unit of 0.8 is inch, about 2cm for the width of one node
    if (gen_max_size <= 16 & label_max_width < 0.8) {
      label_max_width <- 0.8
    }
    if ((label_max_width * gen_max_size) < pdf_max_width) {
      best_cex <- cexs[i] * 0.65
      break
    }
  }

  if (!outline & best_cex == 0 & is.null(cex)) {
    stop("Too many individuals per generation; use compact=TRUE or outline=TRUE.")
  }

  # If user provides cex, recalculate the label_max_width to ensure correct canvas size
  if (!outline && !is.null(cex)) {
    label_max_width <- max(strwidth(max_strwidth_label, cex = cex, units = "inches"), na.rm = TRUE)
    if (gen_max_size <= 16 && label_max_width < 0.8) {
      label_max_width <- 0.8
    }
    if ((label_max_width * gen_max_size) >= pdf_max_width) {
      warning("Provided cex likely exceeds PDF width; labels may clip. Try compact=TRUE or outline=TRUE.")
    }
  }

  #=== Generating the hierarchy layout of all nodes using the sugiyama algorithm =======
  hgap <- round(1 / gen_max_size, 8)
  gen_num <- max(real_node$gen, na.rm = TRUE)
  max_layer <- max(ped_igraph$node$layer, na.rm = TRUE)
  g <- graph_from_data_frame(ped_igraph$edge, directed = TRUE, ped_igraph$node)
  # Map vertex ids to their layer indices directly (preserve reverse layer order).
  layer_idx <- ped_igraph$node[match(V(g)$name, as.character(id)), layer]
  layer_idx <- max_layer - layer_idx + 1
  l <- layout_with_sugiyama(g,
        layers = layer_idx,
        hgap = hgap,
        maxiter = 120,  # Sugiyama layout max iterations
        attributes = "all")$layout
  l <- norm_coords(l,
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1)
  ped_igraph$node <- cbind(ped_igraph$node, x = l[, 1], y = l[, 2])


  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  # Repelling duplicated x positions
  for (i in 1:gen_num) {
    v_rank <- rank(real_node[gen == i, x], na.last = TRUE, ties.method = "first")
    x_sorted <- round(sort(real_node[gen == i, x]), fixed_digits)
    x_new <- repeloverlap(x_sorted)
    real_node[gen == i, x := x_new[v_rank]]
  }
  ped_igraph$node[nodetype %in% c("real", "compact")] <- real_node

  #=== Adjusting space between two nodes (individuals) in x axis for each generation ==
  if ((!outline) & gen_max_size >= 2) {
    x_stats_gen <- real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)), by = gen]
    meanspace <- NULL
    x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
    l_x_distance <- diff(range(l[, 1], na.rm = TRUE))
    max_gen_mean_space <- max(x_stats_gen$meanspace, na.rm = TRUE)
    
    for (i in 1:gen_num) {
      v_rank <- rank(real_node[gen == i, x], na.last = TRUE, ties.method = "first")
      node_num <- length(v_rank)
      if (node_num >= 2) {
        x_distance_1 <- diff(range(real_node[gen == i, x], na.rm = TRUE))
        x_distance_2 <- min((node_num - 1) * max_gen_mean_space, l_x_distance)
        mean_space <- round(x_distance_2 / (node_num - 1), 8)
        x_min <- max(0, min(real_node[gen == i, x], na.rm = TRUE) - (x_distance_2 - x_distance_1))
        x_new <- x_min + seq(from = 0, to = node_num - 1) * mean_space
        real_node[gen == i, x := x_new[v_rank]]
      }
    }
    ped_igraph$node[nodetype %in% c("real", "compact")] <- real_node
  }

  #=== Matching a virtual node's x pos to the smallest position of the full-sib =======
  # A virtual node is a tie between two parents and their progenies
  virtual_node <- ped_igraph$node[nodetype %in% c("virtual")]
  if (nrow(virtual_node) > 0) {
    # Update x via join to preserve row order and column types.
    real_family_min_x <-
      real_node[, .(minx = min(x, na.rm = TRUE)), by = .(gen, familylabel)]
    virtual_node[real_family_min_x, x := i.minx, on = .(gen, familylabel)]
  }
  ped_igraph$node[nodetype %in% c("virtual")] <- virtual_node
  l[, 1] <- ped_igraph$node[, x]


  #=== Rescalling canvas' size, node's size and edge's size ==============================
  # calculate the width of each node: inch
  node_width_s <- label_max_width
  if (!outline) {
    if (gen_max_size >= 2) {
      x_stats_gen <- real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)), by = gen]
      x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
      min_node_space <- min(x_stats_gen[meanspace > 0, meanspace], na.rm = TRUE)
      f <- if (max(x_stats_gen$N) <= 16) 3 * round(node_width_s / min_node_space, 8) else round(node_width_s / min_node_space, 8)
    } else {
      f <- 1
    }
    canvas_width_s <- max(f * l[, 1], na.rm = TRUE) - min(f * l[, 1], na.rm = TRUE) + 6 * node_width_s
  } else {
    # Outline 模式：节点尺寸设为极小值，以显示整体轮廓
    # 逻辑参考原代码：使用 label_max_width 计算画布宽度，但 node_width_s 设为 0.0001
    canvas_width_s <- min(label_max_width * gen_max_size, 200)
    node_width_s <- 0.0001
    
    if (0.0001 * gen_max_size > 200) {
      stop("Outline view unavailable: too many nodes in one generation.")
    }
  }

  # 限制画布宽度在 10-200 英寸之间
  canvas_width_s <- max(10, min(200, canvas_width_s))
  
  # 计算高度：取 (世代高度占用) 和 (宽度 * 0.618) 的较大者，但不超过 200
  canvas_height <- max(gen_num * node_width_s + 3 * node_width_s, canvas_width_s * 0.618)
  
  canvas_height <- min(200, canvas_height)

  node_size <- round(node_width_s * 100 / canvas_width_s, 8)
  V(g)[nodetype %in% c("real", "compact")]$size <- node_size
  
  # Initialize label.cex to avoid NAs or NULL issues
  V(g)$label.cex <- 1
  
  if (outline) {
    node_size <- 0.0001
    V(g)[nodetype %in% c("real", "compact")]$size <- node_size
    # Clear labels for non-highlighted nodes
    V(g)[!V(g)$highlighted]$label <- ""
    
    # Make highlighted nodes visible even in outline mode
    if (any(V(g)$highlighted, na.rm = TRUE)) {
      # Use the original label_max_width to determine a visible size
      highlight_node_size <- round(label_max_width * 100 / canvas_width_s, 8)
      V(g)[V(g)$highlighted & nodetype %in% c("real", "compact")]$size <- highlight_node_size
      
      # Set label cex for highlighted nodes. 
      # If best_cex is 0 (too many nodes), use a default visible cex.
      h_cex <- if (is.null(cex)) {
        if (best_cex > 0) best_cex else 0.6
      } else {
        cex
      }
      V(g)[V(g)$highlighted & nodetype %in% c("real", "compact")]$label.cex <- h_cex
    }
  } else {
    V(g)[nodetype %in% c("real", "compact")]$label.cex <- if (is.null(cex)) best_cex else cex
  }
  
  # Use the maximum node size for edge and margin calculations
  max_node_size <- max(V(g)$size, na.rm = TRUE)

  return(list(
    g = g,
    layout = l,
    canvas_width = canvas_width_s,
    canvas_height = canvas_height,
    node_size = max_node_size,
    best_cex = best_cex
  ))
}

#' Internal plotting function for pedigree igraph
#' @param g An igraph object.
#' @param l A layout matrix.
#' @param node_size Numeric node size.
#' @param ... Additional arguments passed to \code{plot.igraph}.
#' @keywords internal
plot_ped_igraph <- function(g, l, node_size, ...) {
  # Scale attributes directly on the graph to avoid overriding them with global plot arguments.
  # This ensures that edges with arrow.size=0 (like those pointing to virtual nodes) stay arrowless.
  
  # 1. Scale vertex frame width
  v_frame_width <- if (!is.null(V(g)$frame.width)) V(g)$frame.width else 0.2
  # Use a minimum scaling factor to keep borders visible in outline mode
  lwd_scale <- max(if (node_size < 2.5) node_size / 2.5 else 1, 0.1)
  V(g)$frame.width <- v_frame_width * lwd_scale
  
  # 2. Scale edge width
  e_width_base <- if (!is.null(E(g)$width)) E(g)$width else 1
  # Ensure edges are visible even with tiny node_size (e.g. in outline mode)
  E(g)$width <- max(e_width_base * (node_size * 0.15), 0.1)
  
  # 3. Scale arrows
  e_arrow_size_base <- if (!is.null(E(g)$arrow.size)) E(g)$arrow.size else 1
  E(g)$arrow.size <- e_arrow_size_base * (node_size * 0.002)
  
  e_arrow_width_base <- if (!is.null(E(g)$arrow.width)) E(g)$arrow.width else 1
  E(g)$arrow.width <- e_arrow_width_base * (node_size * 0.006)
  
  plot_args <- list(
    x = g,
    rescale = FALSE,
    xlim = c(0 - node_size / 100, 1 + node_size / 100),
    ylim = c(1 + node_size / 100, 0 - node_size / 100),
    layout = l,
    asp = 0
  )
  
  # Override with user-supplied arguments if any
  user_args <- list(...)
  if (length(user_args) > 0) {
    plot_args <- utils::modifyList(plot_args, user_args)
  }
  
  do.call(igraph::plot.igraph, plot_args)
}

ped2igraph <- function(ped, compact = FALSE, highlight = NULL, trace = FALSE, showf = FALSE) {
  if (nrow(ped) == 0) {
    return(list(
      node = data.table(id = integer(), nodetype = character(), gen = integer(), layer = numeric(), label = character()),
      edge = data.table(from = integer(), to = integer())
    ))
  }
  
  # 1. Inject missing parents (for subsetted pedigrees)
  ped_new <- inject_missing_parents(ped)
  
  # 2. Resolve highlight IDs and relatives
  highlight_info <- get_highlight_ids(ped_new, highlight, trace)
  h_ids <- highlight_info$all_ids
  
  # 3. Prepare initial node table
  ped_node <- prepare_initial_nodes(ped_new)
  
  # 4. Compact pedigree (if requested)
  ped_node <- compact_pedigree(ped_node, compact, h_ids)
  
  # 5. Generate edges and virtual nodes
  graph_struct <- generate_graph_structure(ped_node, h_ids)
  ped_node <- graph_struct$node
  ped_edge <- graph_struct$edge
  
  # 6. Apply styles (colors, shapes, highlighting)
  ped_node <- apply_node_styles(ped_node, highlight_info)
  
  # 7. Finalize graph (reindex IDs, set edge colors)
  result <- finalize_graph(ped_node, ped_edge, h_ids, showf)
  
  return(result)
}

# --- Helper Functions ---

inject_missing_parents <- function(ped) {
  ped_new <- copy(ped)
  all_ind_nums <- ped_new$IndNum
  
  # Find missing Sires
  missing_sire_nums <- setdiff(unique(ped_new$SireNum), c(0, NA, all_ind_nums))
  if (length(missing_sire_nums) > 0) {
    sire_info <- unique(ped_new[SireNum %in% missing_sire_nums, .(SireNum, Sire)])
    new_sires <- data.table(
      Ind = sire_info$Sire,
      Sire = NA_character_,
      Dam = NA_character_,
      Sex = "male",
      Gen = min(ped_new$Gen, na.rm=TRUE) - 1,
      IndNum = sire_info$SireNum,
      SireNum = 0L,
      DamNum = 0L
    )
    extra_cols <- setdiff(names(ped_new), names(new_sires))
    if (length(extra_cols) > 0) new_sires[, (extra_cols) := NA]
    ped_new <- rbind(ped_new, new_sires, fill = TRUE)
  }

  # Find missing Dams
  missing_dam_nums <- setdiff(unique(ped_new$DamNum), c(0, NA, all_ind_nums))
  if (length(missing_dam_nums) > 0) {
    dam_info <- unique(ped_new[DamNum %in% missing_dam_nums, .(DamNum, Dam)])
    new_dams <- data.table(
      Ind = dam_info$Dam,
      Sire = NA_character_,
      Dam = NA_character_,
      Sex = "female",
      Gen = min(ped_new$Gen, na.rm=TRUE) - 1,
      IndNum = dam_info$DamNum,
      SireNum = 0L,
      DamNum = 0L
    )
    extra_cols <- setdiff(names(ped_new), names(new_dams))
    if (length(extra_cols) > 0) new_dams[, (extra_cols) := NA]
    ped_new <- rbind(ped_new, new_dams, fill = TRUE)
  }
  
  return(ped_new)
}

get_highlight_ids <- function(ped, highlight, trace) {
  focal_ids <- NULL
  relative_ids <- NULL
  
  # Default colors
  colors <- list(
    focal_frame = "#9c27b0",
    focal_fill = NULL,      # Default: keep original sex-based color
    rel_frame = "#ab47bc",
    rel_fill = NULL         # Default: keep original sex-based color
  )
  
  if (!is.null(highlight)) {
    if (is.list(highlight)) {
      focal_ids <- highlight$ids
      if (!is.null(highlight$frame.color)) colors$focal_frame <- highlight$frame.color
      if (!is.null(highlight$color)) colors$focal_fill <- highlight$color
      if (!is.null(highlight$rel.frame.color)) colors$rel_frame <- highlight$rel.frame.color
      if (!is.null(highlight$rel.color)) colors$rel_fill <- highlight$rel.color
    } else {
      focal_ids <- highlight
    }
    
    # Validate focal_ids
    all_inds <- ped$Ind
    missing_ids <- setdiff(focal_ids, all_inds)
    if (length(missing_ids) > 0) {
      warning("The following highlight IDs were not found in the pedigree: ", paste(missing_ids, collapse = ", "))
    }
    focal_ids <- intersect(focal_ids, all_inds)
    
    # Handle trace
    if (!isFALSE(trace) && length(focal_ids) > 0) {
      trace_dir <- if (isTRUE(trace)) "all" else trace
      relatives_ped <- suppressWarnings(tidyped(ped, cand = focal_ids, trace = trace_dir))
      all_ids <- unique(relatives_ped$Ind)
      relative_ids <- setdiff(all_ids, focal_ids)
    }
  }
  
  list(
    focal = focal_ids,
    relatives = relative_ids,
    all_ids = c(focal_ids, relative_ids),
    colors = colors
  )
}

prepare_initial_nodes <- function(ped) {
  cols <- c("IndNum", "Ind", "SireNum", "DamNum", "Sire", "Dam", "Sex", "Gen")
  if ("Cand" %in% colnames(ped)) cols <- c(cols, "Cand")
  if ("f" %in% colnames(ped)) cols <- c(cols, "f")
  
  ped_node <- ped[, ..cols]
  setnames(ped_node, 
           old = c("IndNum", "Ind", "SireNum", "DamNum", "Sire", "Dam", "Sex", "Gen"),
           new = c("id", "label", "sirenum", "damnum", "sirelabel", "damlabel", "sex", "gen"))
  
  if ("Cand" %in% colnames(ped)) setnames(ped_node, "Cand", "cand")
  
  # Add family info
  max_id <- max(ped_node$id, na.rm = TRUE)
  familylabel = NULL
  ped_node[!(is.na(sirelabel) & is.na(damlabel)),
           familylabel := paste(sirelabel, damlabel, sep = "")]
  
  ped_node[!is.na(familylabel), 
           familynum := .GRP + max_id, 
           by = familylabel]
  ped_node[is.na(familynum), familynum := 0]
  
  nodetype = NULL
  ped_node[, nodetype := "real"]
  
  return(ped_node)
}

compact_pedigree <- function(ped_node, compact, h_ids) {
  if (!compact) return(ped_node)
  
  # Only compact individuals who are not parents of anyone
  sire_dam_label <- unique(c(ped_node$sirelabel, ped_node$damlabel))
  sire_dam_label <- sire_dam_label[!is.na(sire_dam_label)]
  
  ped_node_1 <- ped_node[!(label %in% sire_dam_label) & !is.na(familylabel)]
  
  if (!is.null(h_ids)) {
    ped_node_1 <- ped_node_1[!(label %in% h_ids)]
  }
  
  if (nrow(ped_node_1) > 0) {
    familysize <- NULL
    ped_node_1[, familysize := .N, by = .(familylabel)]
    fullsib_id_DT <- ped_node_1[familysize >= 2]
    
    if (nrow(fullsib_id_DT) > 0) {
      compact_family <- unique(fullsib_id_DT, by = c("familylabel"))
      compact_family[, `:=`(
        label = as.character(familysize),
        nodetype = "compact",
        sex = NA_character_
      )]
      
      next_id <- max(ped_node$id, na.rm = TRUE)
      compact_family[, id := seq.int(next_id + 1, length.out = .N)]
      
      map_dt <- fullsib_id_DT[, .(from_id = id, familylabel)]
      map_dt <- compact_family[map_dt, on = .(familylabel), .(from_id, to_id = id, to_label = label, familylabel)]
      
      ped_node <- ped_node[!(id %in% fullsib_id_DT$id)]
      
      if (nrow(map_dt) > 0) {
        ped_node[map_dt, on = .(sirenum = from_id), `:=`(sirenum = to_id, sirelabel = to_label)]
        ped_node[map_dt, on = .(damnum = from_id), `:=`(damnum = to_id, damlabel = to_label)]
      }
      
      ped_node <- rbind(ped_node, compact_family, fill = TRUE)
      
      # Fix parent references
      compact_lookup <- ped_node[nodetype == "compact", .(compact_id = id, compact_label = label, familylabel)]
      invalid_ids <- setdiff(unique(c(ped_node$sirenum, ped_node$damnum)), ped_node$id)
      invalid_ids <- invalid_ids[invalid_ids != 0 & !is.na(invalid_ids)]
      
      if (length(invalid_ids) > 0 && nrow(compact_lookup) > 0) {
        ped_node[!(sirenum %in% ped_node$id) & !is.na(familylabel), c("sirenum", "sirelabel") := {
          cid <- compact_lookup[.SD, on = .(familylabel), compact_id]
          clab <- compact_lookup[.SD, on = .(familylabel), compact_label]
          list(cid, clab)
        }]
        ped_node[!(damnum %in% ped_node$id) & !is.na(familylabel), c("damnum", "damlabel") := {
          cid <- compact_lookup[.SD, on = .(familylabel), compact_id]
          clab <- compact_lookup[.SD, on = .(familylabel), compact_label]
          list(cid, clab)
        }]
      }
      
      # Recalculate familynum
      max_id <- max(ped_node$id, na.rm = TRUE)
      ped_node[!is.na(familylabel), 
               familynum := .GRP + max_id, 
               by = familylabel]
      ped_node[is.na(familynum), familynum := 0]
    }
  }
  return(ped_node)
}

generate_graph_structure <- function(ped_node, h_ids) {
  # Create edges
  ped_edge <- rbind(
    ped_node[, .(from = id, to = familynum)],
    ped_node[, .(from = familynum, to = sirenum)],
    ped_node[, .(from = familynum, to = damnum)]
  )
  ped_edge <- ped_edge[!(to == 0)]
  ped_edge <- unique(ped_edge)
  ped_edge <- ped_edge[order(from, to)]
  
  width = arrow.size = arrow.width = color = curved = NULL
  edge_default_color <- if (length(h_ids) > 0) "#3333334D" else "#333333"
  
  ped_edge[, ":="(width = 1, arrow.size = 1, arrow.width = 1, arrow.mode = 2, 
                  color = edge_default_color, curved = 0.10)]
  
  # Create virtual nodes
  virtual_nodes <- unique(ped_node[familynum > 0, .(
    id = familynum,
    familylabel,
    label = familylabel,
    sirenum,
    damnum,
    sirelabel,
    damlabel,
    gen,
    familynum,
    highlighted = FALSE
  )])
  
  ped_node <- rbind(ped_node, virtual_nodes, fill = TRUE)
  ped_node[is.na(nodetype), nodetype := "virtual"]
  
  layer = NULL
  ped_node[nodetype %in% c("real", "compact"), layer := 2 * gen - 1]
  ped_node[nodetype %in% c("virtual"), layer := 2 * (gen - 1)]
  
  list(node = ped_node, edge = ped_edge)
}

apply_node_styles <- function(ped_node, highlight_info) {
  shape = frame.color = color = size = label.color = frame.width = label.font = NULL
  
  # Default styles
  ped_node[, ":="(shape = "circle", frame.color = "#7fae59", color = "#9cb383", size = 15, 
                  label.color = "#0d0312", frame.width = 0.2, label.font = 1)]
  ped_node[nodetype == "compact", shape := "square"]
  
  ped_node[sex == "male", ":="(frame.color = "#119ecc", color = "#119ecc")]
  ped_node[sex == "female", ":="(frame.color = "#f4b131", color = "#f4b131")]
  ped_node[is.na(sex) | sex == "unknown", ":="(frame.color = "#9cb383", color = "#9cb383")]
  
  ped_node[nodetype == "virtual", ":="(shape = "none", label = "", size = 0, 
                                       frame.width = 0, color = "#FFFFFF00", frame.color = "#FFFFFF00")]
  
  ped_node[, highlighted := FALSE]
  
  # Apply highlighting
  h_ids <- highlight_info$all_ids
  if (length(h_ids) > 0) {
    focal <- highlight_info$focal
    relatives <- highlight_info$relatives
    colors <- highlight_info$colors
    
    # Mark highlighted nodes
    if (length(relatives) > 0) {
      ped_node[label %in% relatives & nodetype == "real", highlighted := TRUE]
    }
    ped_node[label %in% focal & nodetype == "real", highlighted := TRUE]
    
    h_familynums <- ped_node[highlighted == TRUE, unique(familynum)]
    ped_node[id %in% h_familynums & nodetype == "virtual", highlighted := TRUE]
    
    # Fade non-highlighted nodes
    fade_cols <- function(x) ifelse(nchar(x) == 7, paste0(x, "4D"), x)
    
    ped_node[highlighted == FALSE & nodetype %in% c("real", "compact"), `:=`(
      color = fade_cols(color),
      frame.color = fade_cols(frame.color),
      label.color = fade_cols(label.color)
    )]
    
    # Apply specific colors to highlighted nodes
    if (length(relatives) > 0) {
      ped_node[label %in% relatives & nodetype == "real", 
               ":="(frame.color = colors$rel_frame, 
                    frame.width = 0.5)]
      if (!is.null(colors$rel_fill)) {
        ped_node[label %in% relatives & nodetype == "real", color := colors$rel_fill]
      }
    }
    
    if (length(focal) > 0) {
      ped_node[label %in% focal & nodetype == "real", 
               ":="(frame.color = colors$focal_frame, 
                    frame.width = 1,
                    label.font = 2)]
      if (!is.null(colors$focal_fill)) {
        ped_node[label %in% focal & nodetype == "real", color := colors$focal_fill]
      }
    }
  }
  
  return(ped_node)
}

finalize_graph <- function(ped_node, ped_edge, h_ids, showf) {
  # Reindex IDs
  old_ids <- ped_node$id
  ped_node[, id := seq_len(.N)]
  ped_edge[, from := ped_node$id[match(from, old_ids)]]
  ped_edge[, to := ped_node$id[match(to, old_ids)]]
  
  real_max <- max(ped_node[nodetype %in% c("real", "compact")]$id, na.rm = TRUE)
  
  # Set edge colors
  tonodecolor = i.frame.color = i.highlighted = from_highlighted = NULL
  ped_edge[ped_node, ":="(tonodecolor = i.frame.color, to_highlighted = i.highlighted), on = .(to = id)]
  ped_edge[ped_node, from_highlighted := i.highlighted, on = .(from = id)]
  
  if (length(h_ids) > 0) {
    ped_edge[from_highlighted == TRUE | to_highlighted == TRUE, color := "#333333"]
  }
  
  ped_edge[from > real_max, color := tonodecolor]
  ped_edge[from <= real_max & from_highlighted == TRUE, color := "#333333"]
  
  # Style child-to-family edges
  ped_edge[from <= real_max, ":="(curved = 0, arrow.size = 0, arrow.width = 0, arrow.mode = 0)]
  
  # Sort
  new_names_edge <- c(c("from", "to"), setdiff(colnames(ped_edge), c("from", "to", "tonodecolor")))
  ped_edge <- ped_edge[, ..new_names_edge][order(from, to)]
  
  new_names_node <- c("id", setdiff(colnames(ped_node), "id"))
  ped_node <- ped_node[, ..new_names_node][order(layer, id)]
  
  # Add inbreeding info
  if (showf && "f" %in% colnames(ped_node)) {
    ped_node[nodetype %in% c("real", "compact") & !is.na(f) & f > 0, 
             label := paste0(label, "\n[", round(f, 4), "]")]
  }
  
  list(node = ped_node, edge = ped_edge)
}

`:=` = function(...) NULL

#' Repel overlapping nodes on the x-axis
#'
#' \code{repeloverlap} converts repeated x-axis positions to evenly-spaced continuous positions.
#'
#' When multiple nodes share the same x position, this function spreads them evenly
#' within the available space to their neighbors. For example, given
#' x = c(1.2, 2.1, 2.1, 2.1, 3.2, 4.6, 5.7), where 2.1 appears three times,
#' the function calculates spacing: gap = (3.2 - 2.1) / 3 and transforms to
#' c(1.2, 2.1, 2.1 + gap, 2.1 + 2*gap, 3.2, 4.6, 5.7).
#'
#' @param x A numeric vector of x positions with possible repeated values
#' @return A numeric vector of x positions with unique values
#' @keywords internal
repeloverlap <- function(x) {
  n <- length(x)
  # Early return if no duplicates
  if (n <= 1 || anyDuplicated(x) == 0) {
    return(x)
  }
  
  # Count occurrences using data.table for efficient aggregation
  x_dt <- as.data.table(x)
  x_counts <- x_dt[, .N, by = x]
  
  # Filter to keep only duplicated positions
  dup_info <- x_counts[N > 1]
  if (nrow(dup_info) == 0) {
    return(x)
  }
  
  # Special case: only one unique value (all identical)
  # Cannot spread - return original
  if (nrow(x_counts) == 1) {
    return(x)
  }
  
  # Sort and pre-allocate
  setorder(dup_info, x)
  unique_pos <- sort(x_counts[, x])
  n_unique <- length(unique_pos)
  n_dup <- nrow(dup_info)
  
  # Pre-allocate result vector with maximum possible size
  max_new_vals <- sum(dup_info$N) - n_dup
  result <- numeric(n_unique + max_new_vals)
  result[seq_len(n_unique)] <- unique_pos
  result_idx <- n_unique
  
  # Process each duplicated position
  for (i in seq_len(n_dup)) {
    dup_val <- dup_info$x[i]
    n_copies <- dup_info$N[i]
    idx <- match(dup_val, unique_pos)
    
    # Special case: last two consecutive duplicate positions
    # They share the space between them proportionally
    if (i == n_dup - 1 && n_dup >= 2) {
      next_dup_val <- dup_info$x[i + 1]
      next_idx <- match(next_dup_val, unique_pos)
      
      # Check if they are consecutive in unique_pos
      if (next_idx == idx + 1) {
        # Combine copies from both positions and spread evenly
        next_n_copies <- dup_info$N[i + 1]
        total_copies <- n_copies + next_n_copies - 1
        gap <- (next_dup_val - dup_val) / total_copies
        
        # Generate positions for current duplicates (excluding the first)
        n_new1 <- n_copies - 1
        n_new2 <- next_n_copies - 1
        result[(result_idx + 1):(result_idx + n_new1)] <- dup_val + gap * seq_len(n_new1)
        
        # Generate positions for next duplicates (all of them, excluding last)
        result[(result_idx + n_new1 + 1):(result_idx + n_new1 + n_new2)] <- 
          dup_val + gap * seq(n_copies, total_copies - 1)
        
        result_idx <- result_idx + n_new1 + n_new2
        # Skip the next iteration since we handled it here
        break
      }
    }
    
    # Normal case: spread between current and next position
    n_new <- n_copies - 1
    if (n_new > 0) {
      if (idx < n_unique) {
        next_pos <- unique_pos[idx + 1]
        gap <- (next_pos - dup_val) / n_copies
        # Keep first copy at original position, spread the rest
        result[(result_idx + 1):(result_idx + n_new)] <- dup_val + gap * seq_len(n_new)
      } else if (idx > 1) {
        # Last position: spread backward from previous
        prev_pos <- unique_pos[idx - 1]
        gap <- (dup_val - prev_pos) / n_copies
        # Generate all new positions between prev and current
        result[(result_idx + 1):(result_idx + n_new)] <- prev_pos + gap * seq_len(n_new)
      }
      result_idx <- result_idx + n_new
    }
  }
  
  # Return positions (only the used portion)
  return(result[seq_len(result_idx)])
}
