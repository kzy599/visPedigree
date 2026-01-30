#' Styling and finalizing pedigree graph
#' @import data.table
#' @keywords internal
get_highlight_ids <- function(ped, highlight, trace) {
  focal_ids <- NULL
  relative_ids <- NULL
  
  # Default colors
  colors <- list(
    focal_frame = "#9c27b0",
    focal_fill = NULL,
    rel_frame = "#ab47bc",
    rel_fill = NULL
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
    
    all_inds <- ped$Ind
    focal_ids <- intersect(focal_ids, all_inds)
    
    if (!isFALSE(trace) && length(focal_ids) > 0) {
      trace_dir <- if (isTRUE(trace)) "all" else trace
      relatives_ped <- suppressWarnings(tidyped(ped, cand = focal_ids, trace = trace_dir))
      relative_ids <- setdiff(unique(relatives_ped$Ind), focal_ids)
    }
  }
  
  list(focal = focal_ids, relatives = relative_ids, all_ids = c(focal_ids, relative_ids), colors = colors)
}

#' Apply node styles (color, shape, highlighting)
#' @import data.table
#' @keywords internal
apply_node_styles <- function(ped_node, highlight_info) {
  shape = frame.color = color = size = label.color = frame.width = label.font = highlighted = NULL
  
  # Default styles (set all at once)
  ped_node[, `:=`(
    shape = "circle", 
    frame.color = "#7fae59", 
    color = "#9cb383", 
    size = 15, 
    label.color = "#0d0312", 
    frame.width = 0.2, 
    label.font = 1, 
    highlighted = FALSE
  )]
  
  # Using sub-selections is faster than multiple := lines
  ped_node[nodetype == "compact", shape := "square"]
  ped_node[sex == "male", `:=`(frame.color = "#119ecc", color = "#119ecc")]
  ped_node[sex == "female", `:=`(frame.color = "#f4b131", color = "#f4b131")]
  ped_node[is.na(sex) | sex == "unknown", `:=`(frame.color = "#9cb383", color = "#9cb383")]
  
  # Virtual nodes: circle with tiny size and transparency to fix edge gaps
  ped_node[nodetype == "virtual", `:=`(
    shape = "circle", 
    label = "", 
    size = 0.001, 
    frame.width = 0, 
    color = "#FFFFFF00", 
    frame.color = "#FFFFFF00"
  )]
  
  # Apply highlighting
  h_ids <- highlight_info$all_ids
  if (length(h_ids) > 0) {
    focal <- highlight_info$focal
    relatives <- highlight_info$relatives
    colors <- highlight_info$colors
    
    # Use binary search / keys for matching if many ids
    if (length(relatives) > 0) {
      ped_node[label %in% relatives & nodetype == "real", highlighted := TRUE]
    }
    if (length(focal) > 0) {
      ped_node[label %in% focal & nodetype == "real", highlighted := TRUE]
    }
    
    h_familynums <- ped_node[highlighted == TRUE, unique(familynum)]
    ped_node[id %in% h_familynums & nodetype == "virtual", highlighted := TRUE]
    
    # Fade non-highlighted nodes
    fade_cols <- function(x) ifelse(nchar(x) == 7, paste0(x, "4D"), x)
    # Batch update non-highlighted
    ped_node[highlighted == FALSE & nodetype %in% c("real", "compact"), `:=`(
      color = fade_cols(color), 
      frame.color = fade_cols(frame.color), 
      label.color = fade_cols(label.color)
    )]
    
    # Highlighted styles
    if (length(relatives) > 0) {
      ped_node[label %in% relatives & nodetype == "real", `:=`(
        frame.color = colors$rel_frame, 
        frame.width = 0.5,
        color = if (!is.null(colors$rel_fill)) colors$rel_fill else color
      )]
    }
    
    if (length(focal) > 0) {
      ped_node[label %in% focal & nodetype == "real", `:=`(
        frame.color = colors$focal_frame, 
        frame.width = 1, 
        label.font = 2,
        color = if (!is.null(colors$focal_fill)) colors$focal_fill else color
      )]
    }
  }
  return(ped_node)
}

#' Finalize graph and reindex IDs
#' @import data.table
#' @keywords internal
finalize_graph <- function(ped_node, ped_edge, highlight_info, trace, showf) {
  old_ids <- ped_node$id
  ped_node[, id := seq_len(.N)]
  ped_edge[, from := ped_node$id[match(from, old_ids)]]
  ped_edge[, to := ped_node$id[match(to, old_ids)]]
  
  real_max <- max(ped_node[nodetype %in% c("real", "compact")]$id, na.rm = TRUE)
  
  tonodecolor = i.frame.color = i.highlighted = from_highlighted = NULL
  ped_edge[ped_node, ":="(tonodecolor = i.frame.color, to_highlighted = i.highlighted), on = .(to = id)]
  ped_edge[ped_node, from_highlighted := i.highlighted, on = .(from = id)]
  
  h_ids <- highlight_info$all_ids
  has_trace <- !isFALSE(trace) && length(highlight_info$relatives) > 0
  
  # Default: edges from family nodes to parents follow the parent node color
  ped_edge[from > real_max, color := tonodecolor]
  
  if (length(h_ids) > 0 && has_trace) {
    # When tracing relationships, highlight edges in the path
    # For edges from real nodes to family virtual nodes (individual -> family):
    # Highlight only if the individual (from) is highlighted
    ped_edge[from <= real_max & from_highlighted == TRUE, color := "#333333"]
    
    # For edges from virtual family nodes to parents (family -> parent):
    # Highlight only if BOTH ends are highlighted (family node AND parent)
    ped_edge[from > real_max & from_highlighted == TRUE & to_highlighted == TRUE, color := "#333333"]
  } else if (length(h_ids) == 0) {
    # No highlighting: edges already follow parent node color
  }
  # When highlighting without trace (has_trace == FALSE), keep individual->family edges faded
  
  ped_edge[from <= real_max, ":="(curved = 0, arrow.size = 0, arrow.width = 0, arrow.mode = 0)]
  
  new_names_edge <- c(c("from", "to"), setdiff(colnames(ped_edge), c("from", "to", "tonodecolor")))
  ped_edge <- ped_edge[, ..new_names_edge][order(from, to)]
  
  new_names_node <- c("id", setdiff(colnames(ped_node), "id"))
  ped_node <- ped_node[, ..new_names_node][order(layer, id)]
  
  if (showf && "f" %in% colnames(ped_node)) {
    ped_node[nodetype %in% c("real", "compact") & !is.na(f) & f > 0, 
             label := paste0(label, "\n[", round(f, 4), "]")]
  }
  list(node = ped_node, edge = ped_edge)
}
