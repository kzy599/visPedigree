#' Tidy and prepare a pedigree using graph theory
#'
#' This function takes a pedigree, checks for duplicate and bisexual individuals, detects pedigree loops using graph theory, adds missing founders, assigns generation numbers, sorts the pedigree, and traces the pedigree of specified candidates. If the \code{cand} parameter contains individual IDs, only those individuals and their ancestors or descendants are retained. Tracing direction and the number of generations can be specified using the \code{trace} and \code{tracegen} parameters.
#' 
#' Compared to the legacy version, this function handles cyclic pedigrees more robustly by detecting and reporting loops, and it is generally faster for large pedigrees due to the use of sparse graph algorithms from the \code{igraph} package. Generation assignment can be performed using either a top-down approach (default, aligning founders at the top) or a bottom-up approach (aligning terminal nodes at the bottom).
#'
#' @param ped A data.table or data frame containing the pedigree. The first three columns must be \strong{individual}, \strong{sire}, and \strong{dam} IDs. Additional columns, such as sex or generation, can be included. Column names can be customized, but their order must remain unchanged. Individual IDs should not be coded as "", " ", "0", "*", or "NA"; otherwise, they will be removed. Missing parents should be denoted by "NA", "0", or "*". Spaces and empty strings ("") are also treated as missing parents but are not recommended.
#' @param cand A character vector of individual IDs, or NULL. If provided, only the candidates and their ancestors/descendants are retained.
#' @param trace A character value specifying the tracing direction: "\strong{up}", "\strong{down}", or "\strong{all}". "up" traces ancestors; "down" traces descendants; "all" traces the union of ancestors and descendants. This parameter is only used if \code{cand} is not NULL. Default is "up".
#' @param tracegen An integer specifying the number of generations to trace. This parameter is only used if \code{trace} is not NULL. If NULL or 0, all available generations are traced.
#' @param addgen A logical value indicating whether to generate generation numbers. Default is TRUE, which adds a \strong{Gen} column to the output.
#' @param addnum A logical value indicating whether to generate a numeric pedigree. Default is TRUE, which adds \strong{IndNum}, \strong{SireNum}, and \strong{DamNum} columns to the output.
#' @param inbreed A logical value indicating whether to calculate inbreeding coefficients. Default is FALSE. If TRUE, an \strong{f} column is added to the output.
#' @param genmethod A character value specifying the generation assignment method: "\strong{top}" or "\strong{bottom}". "top" (top-aligned) assigns generations from parents to offspring, starting founders at Gen 1. "bottom" (bottom-aligned) assigns generations from offspring to parents, aligning terminal nodes at the bottom. Default is "top".
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
#' # Trace ancestors of a specific individual within 2 generations
#' tidy_ped_tracegen <- tidyped(simple_ped, cand = "J5X804", trace = "up", tracegen = 2)
#' head(tidy_ped_tracegen)
#' 
#' # Trace both ancestors and descendants for multiple candidates
#' # This is highly optimized and works quickly even on 100k+ individuals
#' cand_list <- c("J5X804", "J3Y620")
#' tidy_ped_all <- tidyped(simple_ped, cand = cand_list, trace = "all")
#' 
#' # Check for loops (will error if loops exist)
#' # tidyped(loop_ped)
#' 
#' # Example with a large pedigree: extract 2 generations of ancestors for 2007 candidates
#' cand_2007 <- big_family_size_ped[Year == 2007, Ind]
#' \donttest{
#' tidy_big <- tidyped(big_family_size_ped, cand = cand_2007, trace = "up", tracegen = 2)
#' summary(tidy_big)
#' }
#'
#' @import data.table
#' @importFrom igraph graph_from_data_frame is_dag topo_sort components subcomponent distances ego degree as_adj_list vcount add_vertices neighbors shortest_paths which_loop ends graph_from_edgelist neighborhood
#' @importFrom stats setNames
#' @importFrom utils head
#' @export
tidyped <- function(ped,
                          cand = NULL,
                          trace = "up",
                          tracegen = NULL,
                          addgen = TRUE,
                          addnum = TRUE,
                          inbreed = FALSE,
                          genmethod = "top",
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

  valid_genmethod <- c("top", "bottom")
  if (length(genmethod) != 1 || !genmethod %in% valid_genmethod) {
    stop(sprintf("The genmethod parameter must be one of: %s", paste(valid_genmethod, collapse = ", ")))
  }
  
  if (!is.null(tracegen)) {
    if (!is.numeric(tracegen) || length(tracegen) != 1 || is.na(tracegen) || !is.finite(tracegen)) {
      stop("The tracegen parameter must be a single numeric value")
    }
    tracegen <- as.integer(tracegen)
  }

  # 2. Data Preparation
  res_prep <- validate_and_prepare_ped(ped)
  ped_dt <- res_prep$ped_dt
  bisexual_parents <- res_prep$bisexual_parents
  
  # 3. Build Initial Graph
  # We need the graph for cycle detection and tracing
  graph_res <- build_ped_graph(ped_dt)
  g <- graph_res$g
  
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
    ped_dt <- assign_ped_generations(g, ped_dt, topo_order, genmethod)
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
  
  # Use match() instead of named vector lookup for speed
  s_from <- match(ped_dt$Sire[!is.na(ped_dt$Sire)], node_names)
  s_to <- which(!is.na(ped_dt$Sire))
  
  d_from <- match(ped_dt$Dam[!is.na(ped_dt$Dam)], node_names)
  d_to <- which(!is.na(ped_dt$Dam))
  
  edge_mat <- matrix(c(s_from, d_from, s_to, d_to), ncol = 2)
  
  # Using graph_from_edgelist is faster usually
  g <- graph_from_edgelist(edge_mat, directed = TRUE)
  
  # Ensure all vertices exist (even isolated ones)
  num_needed <- length(node_names)
  if (vcount(g) < num_needed) {
    g <- add_vertices(g, num_needed - vcount(g))
  }
  
  V(g)$name <- node_names
  
  list(g = g)
}

#' Check for Cycles/Loops in Pedigree
#' @noRd
check_ped_loops <- function(g) {
  if (!is_dag(g)) {
    # 1. Check for self-loops (A -> A)
    # which_loop returns TRUE for edges that are self-loops
    loop_edges <- which(which_loop(g))
    self_loops <- character()
    if (length(loop_edges) > 0) {
      # ends() returns the vertex IDs of edge endpoints
      self_loops <- unique(V(g)[ends(g, loop_edges)[, 1]]$name)
    }
    
    # 2. Check for larger cycles via SCC
    scc <- components(g, mode = "strong")
    cyclic_nodes_idx <- which(scc$csize > 1)
    
    cycle_reports <- character()
    
    if (length(self_loops) > 0) {
      cycle_reports <- c(cycle_reports, paste(self_loops, "->", self_loops, "(self-loop)"))
    }
    
    if (length(cyclic_nodes_idx) > 0) {
      scc_reports <- sapply(cyclic_nodes_idx, function(i) {
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
      cycle_reports <- c(cycle_reports, scc_reports)
    }
    
    # Fallback
    if (length(cycle_reports) == 0) {
      cycle_reports <- "Complex cycle structure detected (try checking strong components)."
    }
    
    stop(paste("Pedigree error! Pedigree loops detected:\n", 
               paste(cycle_reports, collapse = "\n")), call. = FALSE)
  }
}

#' Trace Candidates and Subset Pedigree
#' @noRd
trace_ped_candidates <- function(g, ped_dt, cand, trace, tracegen) {
  # Drop NAs and empty strings from input
  cand_clean <- as.character(cand[!is.na(cand) & cand != "" & cand != " "])
  
  if (length(cand_clean) == 0) {
    stop("The cand parameter contains no valid individual IDs.")
  }

  valid_cand <- cand_clean[cand_clean %in% ped_dt$Ind]
  if (length(valid_cand) == 0) {
    stop("None of the specified candidates were found in the pedigree.")
  }

  if (length(valid_cand) < length(cand_clean)) {
    missing_cands <- setdiff(cand_clean, valid_cand)
    warning(sprintf("The following %d candidates were not found in the pedigree: %s", 
                    length(missing_cands), 
                    paste(head(missing_cands, 5), collapse = ", ")))
  }
  
  # Determine recursion depth
  # order must be numeric and >= 0 for neighborhood()
  order <- if (is.null(tracegen) || tracegen <= 0) vcount(g) else as.integer(tracegen)
  
  # For mode "all", we want the union of ancestors (in) and descendants (out).
  # We should NOT use igraph's mode="all" which is undirected search (connected component).
  if (trace == "all") {
    rel_in <- neighborhood(g, order = order, nodes = valid_cand, mode = "in")
    rel_out <- neighborhood(g, order = order, nodes = valid_cand, mode = "out")
    # Efficiently flatten and get unique IDs from both directions
    keep_indices <- unique(c(
      unlist(lapply(rel_in, as.integer)),
      unlist(lapply(rel_out, as.integer))
    ))
  } else {
    m <- if (trace == "up") "in" else "out"
    reaching_nodes_list <- neighborhood(g, order = order, nodes = valid_cand, mode = m)
    # Efficiently flatten and get unique IDs
    # lapply(list, as.integer) converts vertex sequences to integer IDs
    keep_indices <- unique(unlist(lapply(reaching_nodes_list, as.integer)))
  }
  
  keep_inds <- V(g)$name[keep_indices]
  
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
#' @param genmethod Character, "top" or "bottom".
#' @return A data.table with a \strong{Gen} column added.
#' @noRd
assign_ped_generations <- function(g, ped_dt, topo_order, genmethod) {
  # Get numeric parents for C++ (1-based, 0 for missing)
  # Pedigree is already complete (missing founders added)
  inds <- ped_dt$Ind
  sire_num <- match(ped_dt$Sire, inds, nomatch = 0)
  dam_num <- match(ped_dt$Dam, inds, nomatch = 0)
  
  if (genmethod == "top") {
    gen_vec <- cpp_assign_generations_top(sire_num, dam_num, as.integer(topo_order))
  } else {
    gen_vec <- cpp_assign_generations_bottom(sire_num, dam_num, as.integer(topo_order))
  }
  
  # Isolated nodes check
  in_deg <- degree(g, mode = "in")
  out_deg <- degree(g, mode = "out")
  iso_mask <- (in_deg == 0) & (out_deg == 0)
  if (any(iso_mask)) {
    # Match iso_mask (which is in graph order) to ped_dt order
    # topo_order is a permutation of 1:vcount(g)
    # The C++ gen_vec is in the same order as sire_num/dam_num, which is ped_dt order
    # BUT wait, the graph g was built from ped_dt. So V(g)$name is ped_dt$Ind.
    # So vertex i in the graph corresponds to row i in ped_dt.
    gen_vec[iso_mask] <- 0L
  }
  
  ped_dt[, Gen := gen_vec]
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
