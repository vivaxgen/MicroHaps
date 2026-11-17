#' Plot IBD Network with Dynamic Node and Edge Aesthetics (used to plot Dcifer results)
#'
#' @param data A data.frame containing columns sample_id1, sample_id2, relatedness, and p_value.
#' @param ibd.thres  Numeric. Minimum relatedness value threshold to plot an edge (data less than the minimum threshold will be removed).
#' @param p.val Numeric. Maximum p-value threshold to keep an edge (or data greaterthan the p.value threshold will be removed).
#' @param node.colors.cols Vector of 2 strings. Column names for mapping node colors (NB. the color columns should be 2: 1. for sample_id1 and 2. for sample_id2)).
#' @param node.shape.cols Vector of 2 strings. Column names for mapping node shapes (NB. the shape columns should be 2: 1. for sample_id1 and 2. for sample_id2))..
#' @param edge.thickness Named numeric vector. Key-value pairs for threshold cuts and line widths.
#' @param plot.save.path String. Filepath destination for the output PDF.
#' @param node.shape.name String. Title header text for the shape legend.
#' @param seed Numeric. Random seed calculation state to lock network layout spacing.
#' @param node.size Numeric. Sizing scale parameter passed dynamically to the vertices.
#' @param node.labels Character vector or NULL. Explicit labels for the nodes. If NULL, labels are hidden.
#' @param node.label.cex Numeric. Font character expansion scale (size) for node text labels.
#' @param show.color.legend Logical. If TRUE, renders the bottom-left color classification legend.

#
#' @export
plot_IBD_network <- function(data, 
                             ibd.thres = 0, 
                             p.val = 0.05, 
                             node.colors.cols = NULL, 
                             node.shape.cols = NULL,  
                             edge.thickness = c("0.95" = 2, "0.5" = 1, "0.25" = 0.5), 
                             plot.save.path = "IBD_network_plot.pdf",
                             node.shape.name = "Node Shape",
                             node.size = 8,            
                             show.node.labels = FALSE, # Boolean toggle (Off by default)
                             node.label.cex = 0.8,     # Size controller for labels
                             show.color.legend = TRUE, 
                             seed = 42) {
  
  library(dplyr)
  library(igraph)
  library(RColorBrewer)
  
  # 1. Extract thresholds
  thresh_cuts <- as.numeric(names(edge.thickness))
  sorted_indices <- order(thresh_cuts, decreasing = TRUE)
  thresh_cuts <- thresh_cuts[sorted_indices]
  corresponding_widths <- edge.thickness[sorted_indices]
  
  absolute_min_threshold <- min(thresh_cuts)
  
  # 2. Track Unique Nodes and Filter Edges
  all_nodes <- unique(c(data$sample_id1, data$sample_id2))
  df_edges_filtered <- data %>% 
    filter(p_value < p.val) %>% 
    filter(relatedness >= ibd.thres) %>% 
    filter(relatedness >= absolute_min_threshold)
  
  # 3. Build Graph
  g_net <- graph_from_data_frame(d = df_edges_filtered, vertices = data.frame(name = all_nodes), directed = FALSE)
  node_names <- V(g_net)$name
  
  # 4. Dynamic Shape Mapping
  available_shapes <- c("circle", "square", "rectangle", "csquare", "vrectangle")
  shape_labels <- NULL
  
  if (is.null(node.shape.cols)) {
    V(g_net)$shape <- "circle"
  } else {
    map_shape1 <- data[, c("sample_id1", node.shape.cols[1])]
    colnames(map_shape1) <- c("sample_id", "shape_val")
    map_shape2 <- data[, c("sample_id2", node.shape.cols[2])]
    colnames(map_shape2) <- c("sample_id", "shape_val")
    
    shape_map <- unique(rbind(map_shape1, map_shape2))
    node_shape_vals <- shape_map$shape_val[match(node_names, shape_map$sample_id)]
    shape_labels <- unique(na.omit(node_shape_vals))
    
    shape_indexing <- match(node_shape_vals, shape_labels)
    shape_mapped_indices <- ((shape_indexing - 1) %% length(available_shapes)) + 1
    V(g_net)$shape <- available_shapes[shape_mapped_indices]
  }
  
  # 5. Dynamic Color Mapping
  color_labels <- NULL
  unique_colors_assigned <- NULL
  
  if (is.null(node.colors.cols)) {
    V(g_net)$color <- "skyblue"
  } else {
    map_col1 <- data[, c("sample_id1", node.colors.cols[1])]
    colnames(map_col1) <- c("sample_id", "color_val")
    map_col2 <- data[, c("sample_id2", node.colors.cols[2])]
    colnames(map_col2) <- c("sample_id", "color_val")
    
    color_map <- unique(rbind(map_col1, map_col2))
    node_color_vals <- color_map$color_val[match(node_names, color_map$sample_id)]
    color_labels <- unique(na.omit(node_color_vals))
    num_colors <- length(color_labels)
    
    base_colors <- brewer.pal(min(num_colors, 12), "Set3")
    if (num_colors > 12) {
      base_colors <- colorRampPalette(base_colors)(num_colors)
    }
    
    unique_colors_assigned <- base_colors[1:num_colors]
    V(g_net)$color <- unique_colors_assigned[match(node_color_vals, color_labels)]
  }
  
  # ---------------------------------------------------------------
  # 6. Node Presentation Styles & Automated Text Labels
  # ---------------------------------------------------------------
  V(g_net)$size <- node.size  
  
  if (show.node.labels) {
    V(g_net)$label       <- V(g_net)$name  # Automatically use internal Sample IDs
    V(g_net)$label.cex   <- node.label.cex
    V(g_net)$label.color <- "black"         
    V(g_net)$label.dist  <- 0.75           # Slightly offset text from the shape boundary
  } else {
    V(g_net)$label       <- NA             # Completely hide text
  }
  
  # 7. Edge Styling
  if (ecount(g_net) > 0) {
    edge_shades <- colorRampPalette(c("black", "gray40", "gray70"))(length(thresh_cuts))
    
    E(g_net)$width <- sapply(E(g_net)$relatedness, function(r) {
      match_idx <- which(r >= thresh_cuts)[1]
      if (!is.na(match_idx)) return(corresponding_widths[match_idx])
      return(0.1) 
    })
    
    E(g_net)$color <- sapply(E(g_net)$relatedness, function(r) {
      match_idx <- which(r >= thresh_cuts)[1]
      if (!is.na(match_idx)) return(edge_shades[match_idx])
      return("gray90")
    })
  }
  
  # 8. Render Engine
  pdf(plot.save.path, width = 8, height = 6)
  par(mar = c(1, 1, 1, 1))  
  
  set.seed(seed)
  coords <- layout_with_fr(g_net)
  coords_norm <- coords / max(sqrt(rowSums(coords^2)))
  plot(g_net, layout = coords_norm, margin = 0.2)
  
  usr_coords <- par("usr")   
  box_width <- 0.05 * (usr_coords[2] - usr_coords[1])
  box_height <- 0.06 * (usr_coords[4] - usr_coords[3])
  
  x_left  <- usr_coords[2] - 3.0 * box_width
  x_right <- usr_coords[2] - 2.0 * box_width
  y_top   <- usr_coords[4] - 0.05 * (usr_coords[4] - usr_coords[3])
  
  # 9. Edge Legend
  if (ecount(g_net) > 0) {
    text(x_right + 0.01 * (usr_coords[2] - usr_coords[1]), y_top + 0.5 * box_height, "IBD", adj = 0, cex = 1.2, font = 2)
    y_edge_origin <- y_top - 0.2 * box_height
    for (i in seq_along(thresh_cuts)) {
      y_center <- y_edge_origin - (i - 1) * box_height - box_height / 2
      segments(x0 = (x_left + x_right) / 2, y0 = y_center - box_height / 3, 
               x1 = (x_left + x_right) / 2, y1 = y_center + box_height / 3, 
               col = edge_shades[i], lwd = corresponding_widths[i])
      text(x = x_right + 0.01 * (usr_coords[2] - usr_coords[1]), y = y_center, labels = thresh_cuts[i], adj = 0, cex = 1.1)
    }
  }
  
  # 10. Shape Legend
  if (!is.null(shape_labels)) {
    edge_offset <- if(ecount(g_net) > 0) (length(thresh_cuts) * box_height) else 0
    y_start <- y_top - edge_offset - 0.14 * (usr_coords[4] - usr_coords[3])
    
    text(x_left, y_start + 0.5 * box_height, node.shape.name, adj = 0, cex = 1.2, font = 2)
    y_shape_origin <- y_start - 0.4 * box_height
    for (i in seq_along(shape_labels)) {
      current_shape <- available_shapes[((i - 1) %% length(available_shapes)) + 1]
      y_center <- y_shape_origin - (i-1)*box_height
      if (current_shape == "circle") {
        symbols(x = (x_left + x_right)/2, y = y_center, circles = box_height/3, inches = FALSE, add = TRUE, fg = "black", bg = "gray80")
      } else {
        rect(x_left, y_center - box_height/3, x_right, y_center + box_height/3, border = "black", col = "gray80")
      }
      text(x_right + 0.01 * (usr_coords[2] - usr_coords[1]), y_center, labels = as.character(shape_labels[i]), adj = 0, cex = 1.1)
    }
  }
  
  # 11. Color Legend
  if (show.color.legend && !is.null(color_labels)) {
    legend("bottomleft", legend = color_labels, col = unique_colors_assigned, 
           pch = 16, pt.cex = 1.5, bty = "n", title = "Node Color Groups", 
           title.font = 2, cex = 1.0)
  }
  
  dev.off()
  message(sprintf("Success: Network built. Exactly %d edges drawn.", ecount(g_net)))
}
