#' Visualize a tidy pedigree
#'
#' \code{visped} function draws a graph of a full or compact pedigree.
#'
#' This function takes a pedigree tidied by the \code{\link{tidyped}} function, outputs a hierarchical graph for all individuals in the pedigree. The graph can be shown on the defaulted graphic device and be saved in a pdf file. The graph in the pdf file is a vector drawing, is legible and isn't overlapped. It is especially useful when the number of individuals is big and the width of individual label is long in one generation. This function can draw the graph of a very large pedigree (> 10,000 individuals per generation) by compacting the full-sib individuals. It is very effective for drawing the pedigree of aquatic animal, which usually including many full-sib families per generation in the nucleus breeding population. The outline of a pedigree without individuals' label is still shown if the width of a pedigree graph is longer than the maximum width (200 inches) of the pdf file.
#'
#' In the graph, two shapes and three colors are used. Circle is for individual, square is for family. Dark sky blue means male, dark golden rod means female, dark olive green means unknown sex. For example, one circle with dark sky blue means a male individual; One square with dark golden rod means all female individuals in a full-sib family when \code{compact = TRUE}.
#'
#' @param ped A data.table including the pedigree tidied by the \code{\link{tidyped}} function with the parameter \code{addnum=TRUE}. It is recommended that the pedigree is tidied and pruned by candidates using the \code{\link{tidyped}} function with the non-null parameter \code{cand}.
#' @param compact A logical value indicating whether IDs of full-sib individuals in one generation will be deleted and replaced with the number of full-sib individuals. For example, if there are 100 full-sib individuals in one generation, they will be deleted from the pedigree and be replaced with one individual label of "100" when \code{compact = TRUE}. The default value is FALSE
#' @param outline A logical value indicating whether shapes without label will be shown. A graph of the pedigree without individuals' label is shown when setting \code{outline = TRUE}. It is very useful for viewing the outline of the pedigree and finding the immigrant individuals in each generation when the width of a pedigree graph is longer than the maximum width (200 inches) of the pdf file. The defaulted value is FALSE.
#' @param cex NULL or a numeric value changing the size of individual label shown in the graph. \emph{cex} is an abbreviation of character expansion factor. \code{visped} function will try to guess (\code{cex=NULL}) the matched cex value and returned it in the messages. According to the returned cex of the last run, this parameter should be increased if the label's width is longer than that of the shape in the output pdf file; Contrariwise, this parameter should be decreased if the label's width is shorter than that of the shape in the output pdf file; then rerunning \code{visped} function. The default value is NULL.
#' @param showgraph A logical value indicating whether a plot will be shown in the defaulted graphic device, such as the Plots panel of Rstudio. It is useful for quick viewing of the pedigree graph without opening the pdf file. However, the graph on the defaulted graphic device may be not legible, such as overlapped labels, aliasing lines due to the restricted width and height. It's a good choice to set \code{showgraph = FALSE} when the pedigree is large. The default value is TRUE.
#' @param file NULL or a character value means whether the pedigree graph will be saved in a pdf file. The graph in the pdf file is a legible vector drawing, and labels don't overlap especially when the number of individuals is big and width of the individual label is long in one generation. It is recommended that saving a pedigree graph in the pdf file. The default value is NULL.
#' @return No returned values. The graph will be plotted directly on graphic devices.
#'
#' @examples
#' library(data.table)
#' # Drawing a simple pedigree
#' simple_ped
#' simple_ped_tidy <- tidyped(simple_ped)
#' visped(simple_ped_tidy)
#' # Drawing a simple pedigree of a individual with id of J5X804
#' simple_ped_J5X804_tidy <- tidyped(simple_ped,cand=c("J5X804"))
#' visped(simple_ped_J5X804_tidy)
#' # Drawing the graph in the pdf file
#' visped(simple_ped_J5X804_tidy,file="output.pdf")
#' # Drawing a compact pedigree
#' # The candidates' labels in 2007
#' cand_labels <- big_family_size_ped[(Year == 2007) & (substr(Ind,1,2) == "G8"),Ind]
#' big_ped_tidy <- tidyped(big_family_size_ped,cand=cand_labels)
#' visped(big_ped_tidy,compact=TRUE)
#' visped(big_ped_tidy,compact=TRUE, file="output.pdf")
#' # Individual labels are not shown
#' visped(big_ped_tidy,compact=TRUE, outline=TRUE, file="output.pdf")
#'
#' @import data.table
#' @import igraph
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics strwidth
#' @export
visped <- function(ped,
                   compact = FALSE, outline = FALSE, cex = NULL, showgraph = TRUE, file = NULL) {
  if (is.null(attributes(ped)$tidyped)) {
    stop("The pedigree need to be firstly trimmed by the tidyped() function!")
  } else {
    if (!attr(ped,"tidyped")) {
     stop("The pedigree need to be firstly trimmed by the tidyped() function!")
    }
  }
  ped_new <- copy(ped)

  # Converting pedigree to nodes and edges data.table
  # Reserved digits
  fixed_digits <- 7
  # Digits when calculating
  old_digits <- getOption("digits")
  options(digits = 20)
  # Restore user options even if plotting errors.
  on.exit(options(digits = old_digits), add = TRUE)

  ped_igraph <- ped2igraph(ped_new, compact)
  real_node <- ped_igraph$node[nodetype %in% c("real", "compact")]
  gen_node_num <- real_node[, .N, by = gen]
  gen_max_size <-  max(gen_node_num$N, na.rm = TRUE)


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
      label_max_width = 0.8
    }
    if ((label_max_width * gen_max_size) < pdf_max_width) {
      best_cex <- cexs[i] * 0.7
      break
    }
  }

  if (!outline & best_cex == 0 & is.null(cex)) {
     stop(
       "Too many individuals (>=",
       gen_max_size,
       ") in one generation!!! Two choices:\n", "1. Removing full-sib individuals using the parameter compact = TRUE; or, \n",
       "2. Visualizing all nodes without labels using the parameter outline = TRUE.\n",
       "Rerun visped() function!")
  }
  if (!outline & !is.null(cex)) {
    # Warn when a user-supplied cex likely exceeds PDF width.
    label_max_width_cex <- max(
      strwidth(max_strwidth_label, cex = cex, units = "inches"),
      na.rm = TRUE
    )
    if (gen_max_size <= 16 & label_max_width_cex < 0.8) {
      label_max_width_cex <- 0.8
    }
    if ((label_max_width_cex * gen_max_size) >= pdf_max_width) {
      warning(
        "The provided cex may make the pedigree wider than the PDF maximum width; ",
        "labels or shapes could be clipped. Consider using compact = TRUE or outline = TRUE."
      )
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
  #===Repelling duplicated x positions=================================================
   for (i in 1:gen_num) {
    v_rank <- rank(real_node[gen == i, x], na.last = TRUE, ties.method = "first")
    x_sorted <- sort(real_node[gen == i, x])
    x_sorted <- round(x_sorted, fixed_digits)
    x_new <- repeloverlap(x_sorted)
    real_node[gen == i, x := x_new[v_rank]]
   }
  ped_igraph$node[nodetype %in% c("real", "compact")] <- real_node

  #=== Adjusting space between two nodes (individuals) in x axis for each generation ==
  if ((!outline) & gen_max_size >= 2) {
    x_stats_gen <-
      real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
                by =gen]
    meanspace = NULL # due to NSE notes in R CMD check
    x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
    l_x_range <- range(l[, 1], na.rm = TRUE)
    l_x_distance <- diff(l_x_range)
    max_gen_mean_space <- max(x_stats_gen$meanspace, na.rm = TRUE)
    for (i in 1:gen_num) {
      v_rank <- rank(real_node[gen == i, x], na.last = TRUE, ties.method = "first")
      node_num <- length(v_rank)
      if (node_num >= 2) {
        x_distance_1 <- diff(range(real_node[gen == i, x], na.rm = TRUE))
        x_distance_2 <- (node_num - 1) * max_gen_mean_space
        if (x_distance_2 > l_x_distance) {
          x_distance_2 <- l_x_distance
        }
        mean_space <- round(x_distance_2 / (node_num - 1), 8)
        a <- x_distance_2 - x_distance_1
        x_min <- min(real_node[gen == i, x], na.rm = TRUE)
        b <- x_min - a
        if (b > 0) {
          x_min <- b
        } else {
          x_min <- 0
        }
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
  # l[,1] <- ped_igraph$node[match(V(g)$name,as.character(id)),x]
  l[, 1] <- ped_igraph$node[, x]


  #=== Rescalling canvas' size, node's size and edge's size ==============================
  # calculate the width of each node: inch
  node_width_s <- label_max_width
  if (!outline) {
    if (gen_max_size >= 2) {
      x_stats_gen <-
        real_node[, .(.N, range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
                  by = gen]
      x_stats_gen[range > 0 & N > 1, ":="(meanspace = range / (N - 1))]
      min_node_space <-
        min(x_stats_gen[meanspace > 0, meanspace], na.rm = TRUE)
      f <- 1
      if (max(x_stats_gen$N) <= 16) {
        # increase large space between two nodes when the node number is small
        f <- 3 * round(node_width_s / min_node_space, 8)
      } else {
        # keep small space between two nodes when the node number are big
        f <- round(node_width_s / min_node_space, 8)
      }
    } else {
      f <- 1
    }
    x_f <- f * (l[, 1])
    # Setting the width of the canvas
    # Adding extra 6 node width to the canvas's width to decrease node size
    # when only have one node in one layer
    # because node size is equal to the percentage of node width to the canvas width
    canvas_width_s <- max(x_f, na.rm = TRUE) - min(x_f, na.rm = TRUE) + 6 * node_width_s
  }

  # Finding better node_width when the number of nodes is not too big.
  if (outline) {
    # Maybe < 20,000,000 nodes could be shown in one generation
    node_width_s <- 0.0001
    node_width_v <- seq(from=label_max_width, to=node_width_s, by=-0.0001)
    canvas_width_v <- node_width_v * gen_max_size
    if (min(canvas_width_v,na.rm = TRUE) > pdf_max_width) {
      stop("The outline of the pedigree is not shown due to too many nodes in one generation")
    }
    canvas_width_v <- sort(canvas_width_v)
    canvas_width_s <- max(canvas_width_v[canvas_width_v < pdf_max_width])
  }

  # The maximum width or height for a pdf file is 200 inch
  pdf_maximum_width <- pdf_maximum_height <- pdf_max_width
  if (canvas_width_s > pdf_maximum_width) {
    canvas_width_s <- pdf_maximum_width
  }
  if (canvas_width_s <= 10) {
    canvas_width_s <- 10
  }


  canvas_height <- canvas_width_s * 0.618  # Golden ratio

  # inch
  gen_height <- 0.618  # Golden ratio for generation height
  if (canvas_height < gen_num * (node_width_s) + 3 * node_width_s) {
    canvas_height <- gen_num * (node_width_s) + 3 * node_width_s
  }
  if (canvas_height > pdf_maximum_height) {
    canvas_height <- pdf_maximum_height
  }

  # node_size is a percentage of the width of node to graph
  node_size <- round(node_width_s * 100 / canvas_width_s, 8)
  edge_size <- node_size * 0.001
  edge_arrow_size <- node_size * 0.002
  edge_arrow_width <- node_size * 0.006
  V(g)$size[V(g)$nodetype %in% c("real", "compact")] = node_size
  if (outline) {
    V(g)$label <- ""
  }
  if (!outline) {
    if (is.null(cex) & (best_cex > 0)) {
      V(g)$label.cex[V(g)$nodetype %in% c("real", "compact")] = best_cex
    } else {
      V(g)$label.cex[V(g)$nodetype %in% c("real", "compact")] = cex
    }
    # Avoid NA label.cex values on virtual nodes.
    V(g)$label.cex[is.na(V(g)$label.cex)] <- 1
  }
  E(g)$size = edge_size
  E(g)$arrow.size = edge_arrow_size
  E(g)$arrow.width = edge_arrow_width

  #===Draw the pedigree================================================================
  if (showgraph) {
    plot.igraph(
      g,
      rescale = FALSE,
      xlim = c(0 - node_size / 100, 1 + node_size / 100),
      ylim = c(1 + node_size / 100, 0 - node_size / 100),
      layout = l,
      asp = 0
    )
  }
  if (!is.null(file)) {
    pdf(file = file,
        width = canvas_width_s,
        height = canvas_height)
    # Ensure the device is closed even if plotting fails.
    on.exit(dev.off(), add = TRUE)
    plot.igraph(
      g,
      rescale = FALSE,
      xlim = c(0 - node_size / 100, 1 + node_size / 100),
      ylim = c(1 + node_size / 100, 0 - node_size / 100),
      layout = l,
      asp = 0
    )
    message(paste("The vector drawing of the pedigree is saved in the ",
            getwd(),
            "/",file," file",sep=""))
  }
  if (!outline) {
    if (is.null(cex)) {
      message(paste("The cex for individual label is ", best_cex, ".", sep = ""))
    } else {
      message(paste("The cex for individual label is ", cex, ".", sep = ""))
    }
    message("Please decrease or increase the value of the parameter cex if the label's width is longer or shorter than that of the circle or square in the graph.")

  }
  if (is.null(file)) {
    message("It is recommended that the pedigree graph is saved in the pdf file using the parameter file")
    message("The graph in the pdf file is a vector drawing: shapes, labels and lines are legible; shapes and labels aren't overlapped.")
  }
  invisible(NULL)
}

ped2igraph <- function(ped,compact=TRUE) {
  ped_new <- copy(ped)
  ped_col_names <- colnames(ped_new)
  # There is the Cand column in the pedigree if it is traced by the tidyped function
  if (c("Cand") %in% ped_col_names) {
    ped_node <-
      ped_new[, .(
        id = IndNum,
        label = Ind,
        sirenum = SireNum,
        damnum = DamNum,
        sirelabel = Sire,
        damlabel = Dam,
        cand = Cand,
        sex = Sex,
        gen = Gen
      )]
  } else {
    ped_node <-
      ped_new[, .(
        id = IndNum,
        label = Ind,
        sirenum = SireNum,
        damnum = DamNum,
        sirelabel = Sire,
        damlabel = Dam,
        sex = Sex,
        gen = Gen
      )]
  }

  max_id <- max(ped_node$id,na.rm = TRUE)

  # Adding two new columns family label (column name: familylabel) and it's numeric id
  # (column name: familynum) in the ped_node
  familylabel = NULL # due to NSE notes in R CMD check
  ped_node[!(is.na(sirelabel) &
               is.na(damlabel)),
           familylabel := paste(sirelabel, damlabel, sep = "")]
  family_label <- unique(ped_node$familylabel)
  family_label <- family_label[!is.na(family_label)]
  family_num <-
    setDT(list(
      familynum = seq(
        from = max(ped_node$id,na.rm=TRUE) + 1,
        to = max(ped_node$id,na.rm=TRUE) + length(family_label)
      ),
      familylabel = family_label
    ))
  ped_node <-
    merge(ped_node,
          family_num,
          by = c("familylabel"),
          all.x = TRUE)
  ped_node[is.na(familynum), familynum := 0]

  # There will be three node types in the ped_note, including real, compact, and virtual.
  # Real nodes are all individuals in the pedigree.
  # Compact nodes are full-sib individuals with parents, but without progeny,
  # they exist only when the "compact" parameter is TRUE
  nodetype = NULL # due to NSE notes in R CMD check
  ped_node[,nodetype:="real"]

  #=== Compact the pedigree============================================================
  # Full-sib individuals with parents but without progeny will be deleted from ped_note.
  # count individuals by family and sex as a number of node  replace full-sib individuals
  if (compact) {
    # Finding the individuals with parents, but without progeny
    sire_dam_label <- unique(c(ped_node$sirelabel,ped_node$damlabel))
    sire_dam_label <- sire_dam_label[!is.na(sire_dam_label)]
    ped_node_1 <- ped_node[!(label %in% sire_dam_label)]

    # Moreover, finding full-sib individuals
    familysize <- NULL
    ped_node_1[,familysize:=.N,by=.(familylabel,sex)]
    if (max(ped_node_1$familysize,na.rm=TRUE)>=2) {
      # The full-sib individuals in a family will be compacted if the family size >= 2
      fullsib_id_DT <- ped_node_1[familysize >=2]
      fullsib_ids <- fullsib_id_DT$id
      familylabelsex = NULL # due to NSE notes in R CMD check
      fullsib_id_DT[,familylabelsex:=paste(familylabel,sex,sep="")]
      # Generating a compact family dataset, only including maximum three individuals for
      # each family: male, female and unknown sex individuals
      fullsib_family_label_sex <- unique(fullsib_id_DT$familylabelsex)
      compact_family <- fullsib_id_DT[match(fullsib_family_label_sex,familylabelsex)]
      # The compact families' id are the number of individuals by family and sex.
      compact_family[,":="(label=familysize,nodetype="compact")]
      # Deleting full-sib individuals from families with 2 and more full-sib individuals
      ped_node <- ped_node[!(id %in% fullsib_ids)]
      ped_node <- rbind(ped_node,compact_family,fill=TRUE)
    }
  }

  #=== Add virtual nodes between parents and progrenies================================
  # Add id to familynum and familynum to sirenum and damnum as new virtual edges
  ped_edge <-
    rbind(ped_node[, .(from = id, to = familynum)],
          ped_node[, .(from = familynum, to = sirenum)],
          ped_node[, .(from = familynum, to = damnum)])
  ped_edge <- ped_edge[!(to == 0)]
  # Delete duplicated edges from familynum to parents
  ped_edge <- unique(ped_edge)
  ped_edge <- ped_edge[order(from, to)]
  size = arrow.size = arrow.width = color = curved = NULL # due to NSE notes in R CMD check
  # Set edge colors - currently using dark grey
  ped_edge[,":="(size=1,arrow.size=1,arrow.width=1,color="#333333",curved=0.10)]
  # Add familynum as new virtual nodes
  ped_node <-
    rbind(ped_node, unique(ped_node[familynum > 0, .(
      id = familynum,
      familylabel,
      label = familylabel,
      sirenum,
      damnum,
      sirelabel,
      damlabel,
      gen,
      familynum
    )]), fill = TRUE)
  ped_node[is.na(nodetype),nodetype:="virtual"]
  layer = NULL # due to NSE notes in R CMD check
  ped_node[nodetype %in% c("real","compact"),layer:=2*gen-1]
  ped_node[nodetype %in% c("virtual"),layer:=2*(gen-1)]


  #=== Set default shape, size and color for male and female===========================
  # Setting the default attributes of nodes
  # Notes: size = 15 means the width of a circle node accounts for 15% of the whole width
  # of the graph
  shape = frame.color = color = size = label.color = NULL
  ped_node[, ":="(shape = "circle", frame.color="#7fae59", color="#9cb383",size = 15, label.color="#0d0312")]
  ped_node[nodetype %in% c("compact"), ":="(shape="square")]
  # Setting virtual size of nodes to 0.0001
  ped_node[id > max_id,":="(shape="none",label="",size=0)]
  # Setting male and female's color
  ped_node[sex %in% c("male"), ":="(frame.color="#0e8dbb", color = "#119ecc")]
  ped_node[sex %in% c("female"), ":="(frame.color="#e6a11f", color = "#f4b131")]

  # The edge color is same with the color of the it's "to" node.
  min_familynum <- min(family_num$familynum)
  ped_edge <- merge(ped_edge,
                    ped_node[,.(id,tonodecolor=color)],
                    by.x="to", by.y="id",all.x=TRUE)
  ped_edge[from >= min_familynum,":="(color=tonodecolor)]
  ped_edge[from < min_familynum,":="(curved=0)]


  # Sorting the "from" and "to" columns as the first two columns in the ped_edge
  old_names <- colnames(ped_edge)
  new_names <- c(c("from","to"),old_names[!(old_names %in% c("from","to"))])
  ped_edge <- ped_edge[, ..new_names]
  ped_edge <- ped_edge[order(from,to)]


  # Sorting the "id" column as the first column in the ped_node
  old_names <- colnames(ped_node)
  new_names <- c("id",old_names[!(old_names %in% c("id"))])
  ped_node <- ped_node[, ..new_names]
  ped_node <- ped_node[order(layer,id)]

  return(list(node = ped_node, edge = ped_edge))

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
#' @return A numeric vector of x positions with unique values, preserving order
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
  
  # Return sorted positions (only the used portion)
  return(sort(result[seq_len(result_idx)]))
}
