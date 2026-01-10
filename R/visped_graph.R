#' Convert pedigree to igraph structure
#' @import data.table
#' @keywords internal
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

#' Inject missing parents for subsetted pedigrees
#' @param ped A data.table containing pedigree info.
#' @keywords internal
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

#' Prepare initial node table for igraph conversion
#' @param ped A data.table containing pedigree info.
#' @keywords internal
prepare_initial_nodes <- function(ped) {
  cols <- c("IndNum", "Ind", "SireNum", "DamNum", "Sire", "Dam", "Sex", "Gen")
  if ("Cand" %in% colnames(ped)) cols <- c(cols, "Cand")
  if ("f" %in% colnames(ped)) cols <- c(cols, "f")
  
  ped_node <- ped[, ..cols]
  setnames(ped_node, 
           old = c("IndNum", "Ind", "SireNum", "DamNum", "Sire", "Dam", "Sex", "Gen"),
           new = c("id", "label", "sirenum", "damnum", "sirelabel", "damlabel", "sex", "gen"))
  
  if ("Cand" %in% colnames(ped)) setnames(ped_node, "Cand", "cand")
  
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

#' Compact pedigree by merging full siblings
#' @param ped_node A data.table of nodes.
#' @param compact Logical, whether to compact.
#' @param h_ids Highlighted IDs to exempt from compaction.
#' @keywords internal
compact_pedigree <- function(ped_node, compact, h_ids) {
  if (!compact) return(ped_node)
  
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
      
      compact_lookup <- ped_node[nodetype == "compact", .(compact_id = id, compact_label = label, familylabel)]
      if (length(setdiff(unique(c(ped_node$sirenum, ped_node$damnum)), ped_node$id)) > 0 && nrow(compact_lookup) > 0) {
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
      
      max_id <- max(ped_node$id, na.rm = TRUE)
      ped_node[!is.na(familylabel), 
               familynum := .GRP + max_id, 
               by = familylabel]
      ped_node[is.na(familynum), familynum := 0]
    }
  }
  return(ped_node)
}

#' Generate edges and virtual family nodes
#' @param ped_node A data.table of nodes.
#' @param h_ids Highlighted IDs.
#' @keywords internal
generate_graph_structure <- function(ped_node, h_ids) {
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
  
  virtual_nodes <- unique(ped_node[familynum > 0, .(
    id = familynum,
    familylabel,
    label = familylabel,
    sirenum,
    damnum,
    sirelabel,
    damlabel,
    gen,
    familynum
  )])
  
  ped_node <- rbind(ped_node, virtual_nodes, fill = TRUE)
  ped_node[is.na(nodetype), nodetype := "virtual"]
  
  layer = NULL
  ped_node[nodetype %in% c("real", "compact"), layer := 2 * gen - 1]
  ped_node[nodetype %in% c("virtual"), layer := 2 * (gen - 1)]
  
  list(node = ped_node, edge = ped_edge)
}
