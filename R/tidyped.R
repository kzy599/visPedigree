#' Tidy and prepare a pedigree using graph theory
#'
#' This function takes a pedigree, checks for duplicate and bisexual individuals, detects pedigree loops using graph theory, adds missing founders, assigns generation numbers, sorts the pedigree, and traces the pedigree of specified candidates. If the \code{cand} parameter contains individual IDs, only those individuals and their ancestors or descendants are retained. Tracing direction and the number of generations can be specified using the \code{trace} and \code{tracegen} parameters.
#' 
#' Compared to the legacy version, this function handles cyclic pedigrees more robustly by detecting and reporting loops, and it is generally faster for large pedigrees due to the use of sparse graph algorithms from the \code{igraph} package. Generation assignment is performed using a topological sorting-based algorithm that ensures parents are always placed in a generation strictly above their offspring (TGI algorithm).
#'
#' @param ped A data.table or data frame containing the pedigree. The first three columns must be \strong{individual}, \strong{sire}, and \strong{dam} IDs. Additional columns, such as sex or generation, can be included. Column names can be customized, but their order must remain unchanged. Individual IDs should not be coded as "", " ", "0", "*", or "NA"; otherwise, they will be removed. Missing parents should be denoted by "NA", "0", or "*". Spaces and empty strings ("") are also treated as missing parents but are not recommended.
#' @param cand A character vector of individual IDs, or NULL. If provided, only the candidates and their ancestors/descendants are retained.
#' @param trace A character value specifying the tracing direction: "\strong{up}", "\strong{down}", or "\strong{all}". "up" traces ancestors; "down" traces descendants; "all" traces both simultaneously. This parameter is only used if \code{cand} is not NULL. Default is "up".
#' @param tracegen An integer specifying the number of generations to trace. This parameter is only used if \code{trace} is not NULL. If NULL or 0, all available generations are traced.
#' @param addgen A logical value indicating whether to generate generation numbers. Default is TRUE, which adds a \strong{Gen} column to the output.
#' @param addnum A logical value indicating whether to generate a numeric pedigree. Default is TRUE, which adds \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns to the output.
#' @param inbreed A logical value indicating whether to calculate inbreeding coefficients. Default is FALSE. If TRUE, an \strong{f} column is added to the output.
#' @param ... Additional arguments passed to \code{\link{inbreed}}.
#' 
#' @return A \code{tidyped} object (which inherits from \code{data.table}). Individual, sire, and dam ID columns are renamed to \strong{Ind}, \strong{Sire}, and \strong{Dam}. Missing parents are replaced with \strong{NA}. The \strong{Sex} column contains "male", "female", or NA. The \strong{Cand} column is included if \code{cand} is not NULL. The \strong{Gen} column is included if \code{addgen} is TRUE. The \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns are included if \code{addnum} is TRUE. The \strong{f} column is included if \code{inbreed} is TRUE.
#' 
#' @seealso \code{\link{summary.tidyped}}
#'
#' @examples
#' library(visPedigree)
#' library(data.table)
#'
#' # Tidy a simple pedigree
#' tidy_ped <- tidyped(simple_ped)
#' head(tidy_ped)
#'
#' # Trace ancestors of a specific individual
#' tidy_ped_trace <- tidyped(simple_ped, cand = "J5X804", trace = "up")
#' head(tidy_ped_trace)
#' 
#' # Check for loops (will error if loops exist)
#' # tidyped(loop_ped)
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame is_dag topo_sort components subcomponent distances ego degree as_adj_list vcount add_vertices neighbors shortest_paths
#' @export
tidyped <- function(ped,
                          cand = NULL,
                          trace = "up",
                          tracegen = NULL,
                          addgen = TRUE,
                          addnum = TRUE,
                          inbreed = FALSE,
                          ...) {
  # 1. Parameter Validation
  if (missing(ped) || is.null(ped)) {
    stop("The ped parameter cannot be NULL or missing")
  }
  if ((is.data.frame(ped) || is.matrix(ped)) && NROW(ped) == 0) {
    stop("The ped parameter cannot be empty")
  }
  
  # Boolean checks
  for (arg in c("addgen", "addnum", "inbreed")) {
    val <- get(arg)
    if (!is.logical(val) || length(val) != 1 || is.na(val)) {
      stop(sprintf("The %s parameter only is assigned using TRUE or FALSE", arg))
    }
  }

  valid_trace <- c("up", "down", "all")
  if (length(trace) != 1 || !trace %in% valid_trace) {
    stop(sprintf("The trace parameter must be one of: %s", paste(valid_trace, collapse = ", ")))
  }
  
  if (!is.null(tracegen)) {
    if (!is.numeric(tracegen) || length(tracegen) != 1) {
       stop("The tracegen parameter must be a single numeric value.")
    }
  }

  # 2. Data Preparation
  res_prep <- validate_and_prepare_ped(ped)
  ped_dt <- res_prep$ped_dt
  bisexual_parents <- res_prep$bisexual_parents
  
  # 3. Build Initial Graph
  # We need the graph for cycle detection and tracing
  graph_res <- build_ped_graph(ped_dt)
  g <- graph_res$g
  node_map <- graph_res$node_map # Keep if needed for debugging or advanced mapping
  
  # 4. Cycle Detection
  check_ped_loops(g)
  
  # 5. Trace Candidates (if needed)
  if (!is.null(cand)) {
    ped_dt <- trace_ped_candidates(g, ped_dt, cand, trace, tracegen)
    
    # Rebuild graph for the subsetted pedigree
    # This is necessary because generation assignment depends on the structure of the subset
    graph_res <- build_ped_graph(ped_dt)
    g <- graph_res$g
  }
  
  # 6. Generation Assignment
  # We need topo_order for sorting even if addgen=FALSE
  topo_order <- as.integer(topo_sort(g))
  
  if (addgen) {
    ped_dt <- assign_ped_generations(g, ped_dt, topo_order)
  }
  
  # 7. Sex Inference and Check
  ped_dt <- infer_and_check_sex(ped_dt)
  
  # 8. Sorting and Numeric IDs
  if (addgen) {
    setorder(ped_dt, Gen, Ind)
  } else {
    # If no Gen, preserve topo sort order
    sorted_inds <- V(g)$name[topo_order]
    ped_dt <- ped_dt[match(sorted_inds, Ind)]
  }
  
  if (addnum) {
    ped_dt[, IndNum := .I]
    ped_dt[, SireNum := match(Sire, Ind, nomatch = 0)]
    ped_dt[, DamNum := match(Dam, Ind, nomatch = 0)]
  }
  
  if (!is.null(cand)) {
    ped_dt[, Cand := Ind %in% cand]
  }

  # 9. Final S3 and Inbreeding
  ped_dt <- new_tidyped(ped_dt)
  attr(ped_dt, "bisexual_parents") <- bisexual_parents
  
  if (inbreed) {
    ped_dt <- inbreed(ped_dt, ...)
  }
  
  return(ped_dt)
}

# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Validate and Prepare Pedigree Data
#' @noRd
validate_and_prepare_ped <- function(ped) {
  ped_dt <- if (is.data.table(ped)) copy(ped) else as.data.table(ped)
  
  # Standardize primary column names
  setnames(ped_dt, 1:3, c("Ind", "Sire", "Dam"))
  
  # Ensure character types for IDs
  ped_dt[, c("Ind", "Sire", "Dam") := lapply(.SD, as.character), .SDcols = c("Ind", "Sire", "Dam")]
  
  # Remove records with missing Ind
  if (any(ped_dt$Ind %in% c("", " ", "0", "*", "NA", NA))) {
    warning("Missing values in 'Ind' column; records discarded.")
    ped_dt <- ped_dt[!Ind %in% c("", " ", "0", "*", "NA", NA)]
  }
  
  # Normalize missing parents as NA
  ped_dt[Sire %in% c("", " ", "0", "*", "NA"), Sire := NA_character_]
  ped_dt[Dam %in% c("", " ", "0", "*", "NA"), Dam := NA_character_]
  
  if (all(is.na(ped_dt$Sire)) && all(is.na(ped_dt$Dam))) {
    stop("All parents are missing; no pedigree can be built.")
  }

  # Duplicate checks
  if (anyDuplicated(ped_dt, by = c("Ind", "Sire", "Dam")) > 0) {
    n_duped <- sum(duplicated(ped_dt, by = c("Ind", "Sire", "Dam")))
    warning(sprintf("Removed %d records with duplicate Ind/Sire/Dam IDs.", n_duped))
    ped_dt <- unique(ped_dt, by = c("Ind", "Sire", "Dam"))
  }
  
  if (anyDuplicated(ped_dt, by = "Ind") > 0) {
    stop("Fatal error: Duplicate Ind IDs with different parents found!")
  }

  # Bisexual parent check
  sires <- unique(ped_dt[!is.na(Sire), Sire])
  dams <- unique(ped_dt[!is.na(Dam), Dam])
  bisexual_parents <- sort(intersect(sires, dams))
  
  if (length(bisexual_parents) > 0) {
    warning(paste("Bisexual individuals found (both Sire and Dam):", 
                  paste(bisexual_parents, collapse = ", ")))
  }

  # Add missing founders
  all_parents <- unique(c(sires, dams))
  missing_founders <- setdiff(all_parents, ped_dt$Ind)
  
  if (length(missing_founders) > 0) {
    founder_dt <- data.table(Ind = missing_founders, Sire = NA_character_, Dam = NA_character_)
    ped_dt <- rbind(founder_dt, ped_dt, fill = TRUE)
  }
  
  list(ped_dt = ped_dt, bisexual_parents = bisexual_parents)
}

#' Build igraph Object from Pedigree
#' @noRd
build_ped_graph <- function(ped_dt) {
  node_names <- ped_dt$Ind
  node_map <- setNames(seq_along(node_names), node_names)
  
  # Edges
  s_edges <- ped_dt[!is.na(Sire), .(from = node_map[Sire], to = node_map[Ind])]
  d_edges <- ped_dt[!is.na(Dam), .(from = node_map[Dam], to = node_map[Ind])]
  
  edge_mat <- as.matrix(rbind(s_edges, d_edges))
  
  # Using graph_from_edgelist is faster usually
  g <- graph_from_edgelist(edge_mat, directed = TRUE)
  
  # Ensure all vertices exist (even isolated ones)
  num_needed <- length(node_names)
  if (vcount(g) < num_needed) {
    g <- add_vertices(g, num_needed - vcount(g))
  }
  
  V(g)$name <- node_names
  
  list(g = g, node_map = node_map)
}

#' Check for Cycles/Loops in Pedigree
#' @noRd
check_ped_loops <- function(g) {
  if (!is_dag(g)) {
    scc <- components(g, mode = "strong")
    cyclic_nodes <- which(scc$csize > 1)
    
    cycle_reports <- sapply(cyclic_nodes, function(i) {
      nodes_in_scc <- names(which(scc$membership == i))
      v_start <- nodes_in_scc[1]
      successors <- intersect(names(neighbors(g, v_start, mode = "out")), nodes_in_scc)
      
      if (length(successors) > 0) {
        v_next <- successors[1]
        p <- shortest_paths(g, from = v_next, to = v_start, mode = "out")$vpath[[1]]
        path_names <- c(v_start, names(p))
        paste(path_names, collapse = " -> ")
      } else {
        paste(nodes_in_scc, collapse = ", ")
      }
    })
    
    stop(paste("Pedigree error! Pedigree loops detected:\n", 
               paste(cycle_reports, collapse = "\n")), call. = FALSE)
  }
}

#' Trace Candidates and Subset Pedigree
#' @noRd
trace_ped_candidates <- function(g, ped_dt, cand, trace, tracegen) {
  valid_cand <- cand[cand %in% ped_dt$Ind]
  if (length(valid_cand) == 0) {
    stop("None of the specified candidates were found in the pedigree.")
  }
  
  keep_inds <- character(0)
  
  if (trace == "up") {
    for (c in valid_cand) {
      anc <- subcomponent(g, c, mode = "in")
      if (!is.null(tracegen) && tracegen > 0) {
        dists <- distances(g, v = c, to = anc, mode = "in")
        anc <- anc[as.vector(dists) <= tracegen]
      }
      keep_inds <- union(keep_inds, names(anc))
    }
  } else if (trace == "down") {
    for (c in valid_cand) {
      des <- subcomponent(g, c, mode = "out")
      if (!is.null(tracegen) && tracegen > 0) {
        dists <- distances(g, v = c, to = des, mode = "out")
        des <- des[as.vector(dists) <= tracegen]
      }
      keep_inds <- union(keep_inds, names(des))
    }
  } else if (trace == "all") {
    for (c in valid_cand) {
      anc <- subcomponent(g, c, mode = "in")
      des <- subcomponent(g, c, mode = "out")
      if (!is.null(tracegen) && tracegen > 0) {
        anc <- anc[as.vector(distances(g, v = c, to = anc, mode = "in")) <= tracegen]
        des <- des[as.vector(distances(g, v = c, to = des, mode = "out")) <= tracegen]
      }
      keep_inds <- union(keep_inds, union(names(anc), names(des)))
    }
  }
  
  ped_dt <- ped_dt[Ind %in% keep_inds]
  # Update parents not in set after tracing
  ped_dt[!(Sire %in% keep_inds), Sire := NA_character_]
  ped_dt[!(Dam %in% keep_inds), Dam := NA_character_]
  
  return(ped_dt)
}

#' Assigns individual generation numbers based on topological sorting and parentage.
#'
#' @param g An igraph object representing the pedigree.
#' @param ped_dt A data.table containing the pedigree.
#' @param topo_order Integer vector specifying the topological order of vertices.
#' @return A data.table with a \strong{Gen} column added.
#' @noRd
assign_ped_generations <- function(g, ped_dt, topo_order) {
  # 1. Height Calculation (Bottom-Up Logic)
  
  # Optimization: Skip nodes with no parents (In-Degree 0)
  in_deg <- degree(g, mode = "in")
  has_parents_mask <- in_deg[topo_order] > 0
  
  rev_topo_proc <- rev(topo_order[has_parents_mask])
  
  height_vec <- rep(0L, vcount(g))
  adj_list_in <- as_adj_list(g, mode = "in")
  
  for (v_idx in rev_topo_proc) {
    parents <- adj_list_in[[v_idx]]
    p_indices <- as.integer(parents)
    prior_h <- height_vec[p_indices]
    new_h <- height_vec[v_idx] + 1L
    height_vec[p_indices] <- pmax(prior_h, new_h)
  }
  
  # 2. Align Full Siblings
  # Identify individuals who have at least one parent
  h_dt <- data.table(
    Ind = V(g)$name,
    Height = height_vec,
    Sire = ped_dt$Sire[match(V(g)$name, ped_dt$Ind)],
    Dam = ped_dt$Dam[match(V(g)$name, ped_dt$Ind)]
  )
  
  non_orphans <- h_dt[!is.na(Sire) | !is.na(Dam)]
  
  if (nrow(non_orphans) > 0) {
    non_orphans[, MaxH := max(Height), by = .(Sire, Dam)]
    update_inds <- non_orphans[MaxH > Height]
    if (nrow(update_inds) > 0) {
      update_indices <- match(update_inds$Ind, V(g)$name)
      height_vec[update_indices] <- update_inds$MaxH
    }
  }
  
  gen_map_height <- setNames(height_vec, V(g)$name)
  founders <- names(which(degree(g, mode = "in") == 0))
  
  if (length(founders) > 0) {
    pairs_dt <- unique(ped_dt[!is.na(Sire) & !is.na(Dam), .(Sire, Dam)])
    
    if (nrow(pairs_dt) > 0) {
      pairs_dt[, SireH := gen_map_height[Sire]]
      pairs_dt[, DamH := gen_map_height[Dam]]
      
      updates <- list()
      
      fs_mask <- pairs_dt$Sire %in% founders & pairs_dt$SireH < pairs_dt$DamH
      if (any(fs_mask)) {
        updates$sires <- pairs_dt[fs_mask, .(MaxMateH = max(DamH)), by = Sire]
      }
      
      fd_mask <- pairs_dt$Dam %in% founders & pairs_dt$DamH < pairs_dt$SireH
      if (any(fd_mask)) {
        updates$dams <- pairs_dt[fd_mask, .(MaxMateH = max(SireH)), by = Dam]
      }
      
      if (!is.null(updates$sires)) {
        gen_map_height[updates$sires$Sire] <- updates$sires$MaxMateH
      }
      if (!is.null(updates$dams)) {
        gen_map_height[updates$dams$Dam] <- updates$dams$MaxMateH
      }
      
      height_vec <- gen_map_height
    }
  }
  
  # 3. Propagate Height Down (Parent -> Child)
  # Enforce Child H >= min(Parent H) - 1
  topo_proc <- topo_order[has_parents_mask]
  
  for (v_idx in topo_proc) {
    parents <- adj_list_in[[v_idx]]
    p_indices <- as.integer(parents)
    p_heights <- height_vec[p_indices]
    
    min_p_h <- min(p_heights)
    if (min_p_h - 1L > height_vec[v_idx]) {
      height_vec[v_idx] <- min_p_h - 1L
    }
  }
  
  gen_map_height <- setNames(height_vec, V(g)$name)
  max_h <- if (length(gen_map_height)>0) max(gen_map_height) else 0
  final_gens <- max_h - gen_map_height + 1L
  
  # Identify isolated vertices (no parents, no offspring) and assign Gen 0
  # This matches legacy behavior and allows visped to filter them out
  iso_nodes <- names(which(degree(g, mode = "all") == 0))
  if (length(iso_nodes) > 0) {
    final_gens[iso_nodes] <- 0L
  }
  
  ped_dt[, Gen := final_gens[match(Ind, names(final_gens))]]
  
  return(ped_dt)
}

#' Infer and Check Sex of Individuals
#' @noRd
infer_and_check_sex <- function(ped_dt) {
  if (!("Sex" %in% names(ped_dt))) {
    ped_dt[, Sex := NA_character_]
  } else {
    ped_dt[, Sex := tolower(as.character(Sex))]
    ped_dt[Sex %in% c("", " ", "na"), Sex := NA_character_]
  }
  
  sires <- unique(ped_dt[!is.na(Sire), Sire])
  dams <- unique(ped_dt[!is.na(Dam), Dam])
  
  # Sex conflict check
  sex_conflicts <- ped_dt[(Ind %in% sires & Sex == "female") | (Ind %in% dams & Sex == "male"), Ind]
  if (length(sex_conflicts) > 0) {
    warning("Potential sex conflicts detected for individuals: ", paste(sex_conflicts, collapse = ", "))
  }
  
  # Infer sex from roles
  ped_dt[is.na(Sex) & (Ind %in% sires), Sex := "male"]
  ped_dt[is.na(Sex) & (Ind %in% dams), Sex := "female"]
  
  return(ped_dt)
}
