#' Render pedigree graph using Two-Pass strategy
#' @importFrom igraph V E plot.igraph
#' @importFrom utils modifyList
#' @keywords internal
plot_ped_igraph <- function(g, l, node_size, ...) {
  # ============================================================================
  # SCALING LOGIC FOR RENDERING
  # ============================================================================
  # In 'outline' mode, node_size can be extremely small (e.g., 0.0001).
  # We use a 'scaling_ref' to ensure edges and arrows remain visible.
  # For normal pedigrees, node_size is typically 2-10.
  scaling_ref <- max(node_size, 1.5)
  
  # Scale frame width
  v_frame_width <- if (!is.null(V(g)$frame.width)) V(g)$frame.width else 0.2
  # frame.width should be smaller for tiny nodes but not vanish
  lwd_scale <- if (scaling_ref < 2.5) scaling_ref / 2.5 else 1
  V(g)$frame.width <- v_frame_width * lwd_scale
  
  # Scale edge width
  e_width_base <- if (!is.null(E(g)$width)) E(g)$width else 1
  E(g)$width <- pmax(e_width_base * (scaling_ref * 0.15), 0.001)
  
  # Scale arrows
  e_arrow_size_base <- if (!is.null(E(g)$arrow.size)) E(g)$arrow.size else 1
  calc_arrow_size <- e_arrow_size_base * (scaling_ref * 0.005)
  
  if (median(calc_arrow_size, na.rm = TRUE) < 0.005) { # Lower threshold for tiny plots
    E(g)$arrow.size <- 0
    E(g)$arrow.width <- 0
    E(g)$arrow.mode <- 0
  } else {
    E(g)$arrow.size <- calc_arrow_size
    e_arrow_width_base <- if (!is.null(E(g)$arrow.width)) E(g)$arrow.width else 1
    E(g)$arrow.width <- e_arrow_width_base * (scaling_ref * 0.015)
    if (is.null(E(g)$arrow.mode)) E(g)$arrow.mode <- 2
  }
  
  user_args <- list(...)
  # Margin should be at least large enough to show a bit of the grid
  margin <- max(node_size / 100, 0.02)
  
  # PASS 1: Draw EDGES ONLY
  g_edges <- g
  original_sizes <- V(g_edges)$size
  V(g_edges)$size <- 0.001
  V(g_edges)$color <- NA
  V(g_edges)$frame.color <- NA
  V(g_edges)$frame.width <- 0
  V(g_edges)$label <- NA
  
  plot_args_edges <- list(
    x = g_edges, rescale = FALSE,
    xlim = c(-margin, 1 + margin), ylim = c(1 + margin, -margin),
    layout = l, asp = 0, add = FALSE, vertex.label = NA
  )
  
  if (length(user_args) > 0) {
    safe_args <- user_args[!grepl("^vertex\\.", names(user_args))]
    plot_args_edges <- utils::modifyList(plot_args_edges, safe_args)
    plot_args_edges$vertex.label <- NA
    plot_args_edges$vertex.size <- 0.001
  }
  suppressWarnings(do.call(igraph::plot.igraph, plot_args_edges))
  
  # PASS 2: Draw NODES AND LABELS ONLY
  g_nodes <- g
  V(g_nodes)$size <- original_sizes
  E(g_nodes)$color <- "#FFFFFF00"
  E(g_nodes)$width <- 0
  E(g_nodes)$arrow.size <- 0
  E(g_nodes)$arrow.mode <- 0
  
  plot_args_nodes <- list(
    x = g_nodes, rescale = FALSE,
    xlim = c(-margin, 1 + margin), ylim = c(1 + margin, -margin),
    layout = l, asp = 0, add = TRUE
  )
  if (length(user_args) > 0) plot_args_nodes <- utils::modifyList(plot_args_nodes, user_args)
  suppressWarnings(do.call(igraph::plot.igraph, plot_args_nodes))
}
