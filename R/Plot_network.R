#' @title Plot ceRNA network
#' @description This function plots a ceRNA network, including circRNAs, lncRNAs, miRNAs and mRNAs
#' @param net A data frame that contains information about the nodes in the network. It should have two columns: "From" and "To". Use "UPnetwork", "DOWNnetwork" or "mynetwork". 
#' @param nodeshape The shape of the nodes. Default is "circle". If you want to draw a 3D plot, use "sphere" instead.
#' @param circ_color a character specifying the color for circular RNA nodes
#' @param lnc_color a character specifying the color for long non-coding RNA nodes
#' @param mir_color a character specifying the color for microRNA nodes
#' @param Up_mrna_color a character specifying the color for upregulated mRNA nodes
#' @param Down_mrna_color a character specifying the color for downregulated mRNA nodes
#' @param circ_alpha a numeric value specifying the alpha transparency for circular RNA nodes
#' @param lnc_alpha a numeric value specifying the alpha transparency for long non-coding RNA nodes
#' @param mir_alpha a numeric value specifying the alpha transparency for microRNA nodes
#' @param Up_alpha a numeric value specifying the alpha transparency for upregulated mRNA nodes
#' @param Down_alpha a numeric value specifying the alpha transparency for downregulated mRNA nodes
#' @param circ_size a numeric value specifying the size for circular RNA nodes
#' @param lnc_size a numeric value specifying the size for long non-coding RNA nodes
#' @param mir_size a numeric value specifying the size for microRNA nodes
#' @param mrna_size a numeric value specifying the size for mRNA nodes
#' @param circ_labelsize a numeric value specifying the size for labels of circular RNA nodes
#' @param lnc_labelsize a numeric value specifying the size for labels of long non-coding RNA nodes
#' @param mir_labelsize a numeric value specifying the size for labels of microRNA nodes
#' @param mrna_labelsize a numeric value specifying the size for labels of mRNA nodes
#' @param label_color The color of the node labels. Default is "black".
#' @param edge_color The color of the edges. Default is "grey".
#' @param edge_alpha The alpha value of the edges. Default is 0.75.
#' @param netlayout The layout of the network.  Default is NULL, and the function will calculate the layout using the Kamada-Kawai algorithm.
#' @param pdf.name The name of the PDF file to be saved. Default is "TFnet.pdf".
#' @param pdf.height The height of the PDF file. Default is 15.
#' @param pdf.width The width of the PDF file. Default is 15.
#' @param edge.arrow.size The size of the arrows on the edges. Default is 0.4.
#' @param vertex.label.family The font family of the node labels. Default is "sans".
#' @return A PDF file containing the network plot, and the network object stored in the global environment
#' @export
#' @import igraph
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import circlize
plot_CEnetwork <- function(net,
                           nodeshape = "circle",
                           circ_color = "#459943",
                           lnc_color = "#F0E442",
                           mir_color = "#CC79A7",
                           Up_mrna_color = "#fb8072",
                           Down_mrna_color = "#80b1d3",
                           circ_alpha = 0.5,
                           lnc_alpha = 0.5,
                           mir_alpha = 0.4,
                           Up_alpha = 0.3,
                           Down_alpha = 0.3,
                           circ_size = 10,
                           lnc_size = 10,
                           mir_size = 8,
                           mrna_size = 8,
                           circ_labelsize = 0.7,
                           lnc_labelsize = 0.7,
                           mir_labelsize = 0.7,
                           mrna_labelsize = 0.7,
                           label_color = "black",
                           edge_color = "grey",
                           edge_alpha = 0.75,
                           netlayout = NULL,
                           pdf.name = "TFnet.pdf",
                           pdf.height = 15,
                           pdf.width = 15,
                           edge.arrow.size = 0.4,
                           vertex.label.family = "sans") {
  #  options(device = "RStudio")
  library(igraph)
  all_nodes <- unique(c(net$From, net$To))
  
  nodecolor <- ifelse(all_nodes %in% DEcirc$degs,
                      adjustcolor(circ_color, alpha.f = circ_alpha),
                      ifelse(all_nodes %in% DElnc$degs,
                             adjustcolor(lnc_color, alpha.f = lnc_alpha),
                             ifelse(all_nodes %in% DEmir$degs,
                                    adjustcolor(mir_color, alpha.f = mir_alpha),
                                    ifelse(all_nodes %in% DEm$up,
                                           adjustcolor(Up_mrna_color, alpha.f = Up_alpha),
                                           adjustcolor(Down_mrna_color, alpha.f = Down_alpha)))))
  
  
  nodesize <- ifelse(all_nodes %in% DEcirc$degs,
                     scale(table(DEcirc$degs), center = FALSE) * circ_size + 0.6 ,
                     ifelse(all_nodes %in% DElnc$degs,
                            scale(table(DElnc$degs), center = FALSE) * lnc_size + 0.6,
                            ifelse(all_nodes %in% DEmir$degs,
                                   scale(table(DEmir$degs), center = FALSE) * mir_size + 0.6,
                                   mrna_size)))
  
  labelsize <- ifelse(all_nodes %in% DEcirc$degs,
                      scale(table(DEcirc$degs), center = FALSE) * circ_labelsize + 0.6,
                      ifelse(all_nodes %in% DElnc$degs,
                             scale(table(DElnc$degs), center = FALSE) * lnc_labelsize + 0.6,
                             ifelse(all_nodes %in% DEmir$degs,
                                    scale(table(DEmir$degs), center = FALSE) * mir_labelsize + 0.2,
                                    mrna_labelsize)))
  
  edgecolor <- adjustcolor(edge_color, alpha.f = edge_alpha)
  
  node <- data.frame(nodes = all_nodes,
                     nodesize =  nodesize,
                     nodecolor = nodecolor)
  
  
  edgecolor <- rep(edgecolor, nrow(net))
  
  #定义network
  link <- net
  colnames(link) <- c("From", "To")
  link$color <- edgecolor
  
  network  <- graph_from_data_frame(d=link, vertices=node, directed=T)
  
  V(network)$color <- V(network)$nodecolor
  V(network)$size <- V(network)$nodesize 
  V(network)$label.cex <- labelsize
  V(network)$shape <- nodeshape
  V(network)$label.color <- label_color 
  
  edge.start <- ends(network, es = E(network), names = FALSE)[, 1] 
  edge.col <- V(network)$color[edge.start]
  E(network)$color <- E(network)$color
  
  #设置网络图的布局
  if (is.null(netlayout)) {
    netlayout <- layout_(network, 
                         with_kk())
  }
  
  #将网络对象存储在内存中
  assign("network", network, envir = .GlobalEnv)
  
  
  #保存为 PDF 文件
  pdf(file = pdf.name, height = pdf.height, width = pdf.width)
  plot(network, layout = netlayout, edge.arrow.size = edge.arrow.size,
       vertex.label.family = vertex.label.family, res = c(1200, 1200),
       rescale = TRUE)
  dev.off()
  
}


get_downstream <- function(network, highlight_node) {
  downstream_nodes <- names(neighbors(network, highlight_node, mode = "out"))
  second_downstream_nodes <- lapply(downstream_nodes, function(x) names(neighbors(network, x, mode = "out")))
  second_downstream_nodes2 <- unlist(second_downstream_nodes)
  all_downstream_nodes <- c(highlight_node, downstream_nodes, second_downstream_nodes2)
  return(all_downstream_nodes)
}


get_upstream <- function(network, highlight_node) {
  upstream_nodes <- names(neighbors(network, highlight_node, mode = "in"))
  second_upstream_nodes <- lapply(upstream_nodes, function(x) names(neighbors(network, x, mode = "in")))
  second_upstream_nodes2 <- unlist(second_upstream_nodes)
  all_upstream_nodes <- c(highlight_node, upstream_nodes, second_upstream_nodes2)
  return(all_upstream_nodes)
}

#' @title Plot a ceRNA network with highlighted RNAs of interest
#' @description This function plots a highlighted network with a specified node highlighted in a different color. The function takes as input the network and the node to be highlighted, and outputs a PDF file of the network plot with the highlighted node and its connected edges colored differently.
#' @param highlight_node Your RNA of interest
#' @param network A graph object representing the ceRNA network, created using "plot_CEnetwork" before.
#' @param highlight_circ_color Color of highlighted circRNAs. Default is "#459943"
#' @param highlight_lnc_color Color of highlighted lncRNAs. Default is "#F0E442"
#' @param highlight_mir_color Color of highlighted miRNAs, Default is "#CC79A7"
#' @param highlight_edge_color Color of highlighted edges. Default is "#e41a1c"
#' @param highlight_circ_alpha a numeric value specifying the alpha transparency for highlighted circular RNA nodes
#' @param highlight_lnc_alpha a numeric value specifying the alpha transparency for highlighted long non-coding RNA nodes
#' @param highlight_mir_alpha a numeric value specifying the alpha transparency for highlighted microRNA nodes
#' @param highlight_m_alpha a numeric value specifying the alpha transparency for highlighted  mRNA nodes
#' @param netlayout A layout object representing the layout of the network. Default is NULL, and the function will calculate the layout using the Kamada-Kawai algorithm.
#' @param pdf_height An integer representing the height of the PDF output file in inches. Default is 15.
#' @param pdf_width An integer representing the width of the PDF output file in inches. Default is 15.
#' @param vertex_label_family A character string representing the font family for the vertex labels in the network plot. Default is "sans".
#' @return This function outputs a PDF file of the network plot with the highlighted node and its connected edges colored differently.
#' @export
#' @import igraph
#' @import magrittr
#' @import dplyr
#' @import ggplot2
#' @import circlize
plot_highlighted_CEnetwork <- function(network, 
                                       highlight_node, 
                                       highlight_circ_color = "#459943", 
                                       highlight_lnc_color = "#F0E442", 
                                       highlight_mir_color = "#CC79A7", 
                                       highlight_edge_color = "#e41a1c",
                                       highlight_circ_alpha = 0.9,
                                       highlight_lnc_alpha = 0.9,
                                       highlight_mir_alpha = 0.9,
                                       highlight_m_alpha = 0.9,
                                       netlayout = NULL,
                                       pdf_height = 15, 
                                       pdf_width = 15, 
                                       vertex_label_family = "sans") {
  # 判断 highlight_node 是lnc circ mir OR mRNA.
  if (highlight_node %in% c(DElnc$degs, DEcirc$degs)) {
    # 如果是，将高亮节点设置为 TF 颜色，与之相连的节点设置为目标颜色
    
    all_downstream_nodes <- get_downstream(network, highlight_node)
    
    for (node in all_downstream_nodes) {
      if (node %in% DEmir$degs) {
        V(network)[node]$color <- adjustcolor(highlight_mir_color, alpha.f = highlight_mir_alpha)
      } else if (node %in% DElnc$degs) {
        V(network)[node]$color <- adjustcolor(highlight_lnc_color, alpha.f = highlight_lnc_alpha)
      } else if (node %in% DEcirc$degs) {
        V(network)[node]$color <- adjustcolor(highlight_circ_color, alpha.f = highlight_circ_alpha)
      } else if (node %in% DEm$up) {
        V(network)[node]$color <- adjustcolor("#fb8072", alpha.f = highlight_m_alpha)
      } else {
        V(network)[node]$color <- adjustcolor("#80b1d3", alpha.f = highlight_m_alpha)
      }
    }
    
    
  } else if (highlight_node %in% DEmir$degs) {
    # 如果不是，将高亮节点设置为目标颜色，与之相连的节点设置为 TF 颜色
    connected_nodes <- neighbors(network, highlight_node, mode = "all")
    connected_nodes2 <- neighbors(network, highlight_node, mode = "in")
    connected_nodes3 <- c(highlight_node, names(connected_nodes2))
    connected_nodes4 <- neighbors(network, highlight_node, mode = "out")
    connected_nodes5 <- c(highlight_node, names(connected_nodes4))   
    V(network)[highlight_node]$color <- (adjustcolor(highlight_mir_color, alpha.f = 0.9))
    
    for (node in names(connected_nodes)) {
      if (node %in% DElnc$degs) {
        V(network)[node]$color <- adjustcolor(highlight_lnc_color, alpha.f = highlight_lnc_alpha)
      } else if (node %in% DEcirc$degs) {
        V(network)[node]$color <- adjustcolor(highlight_circ_color, alpha.f = highlight_circ_alpha)
      } else if (node %in% DEm$up) {
        V(network)[node]$color <- adjustcolor("#fb8072", alpha.f = highlight_m_alpha)
      } else {
        V(network)[node]$color <- adjustcolor("#80b1d3", alpha.f = highlight_m_alpha)
      }
    }
    
  } else {
    # 如果不是，将高亮节点设置为目标颜色，与之相连的节点设置为 TF 颜色
    all_upstream_nodes <- get_upstream(network, highlight_node)
    
    for (node in all_upstream_nodes) {
      if (node %in% DEmir$degs) {
        V(network)[node]$color <- adjustcolor(highlight_mir_color, alpha.f = highlight_mir_alpha)
      } else if (node %in% DElnc$degs) {
        V(network)[node]$color <- adjustcolor(highlight_lnc_color, alpha.f = highlight_lnc_alpha)
      } else if (node %in% DEcirc$degs) {
        V(network)[node]$color <- adjustcolor(highlight_circ_color, alpha.f = highlight_circ_alpha)
      } else if (node %in% DEm$up) {
        V(network)[node]$color <- adjustcolor("#fb8072", alpha.f = highlight_m_alpha)
      } else {
        V(network)[node]$color <- adjustcolor("#80b1d3", alpha.f = highlight_m_alpha)
      }
    }
    
  }
  
  # 获取所有与 highlight_node 相连的边
  if (highlight_node %in% c(DElnc$degs, DEcirc$degs)) {
    connected_edges <- incident_edges(network, all_downstream_nodes, mode = "out")
    E(network)[unlist(connected_edges)]$color <- highlight_edge_color
  } else if (highlight_node %in% DEmir$degs) {
    connected_edges <- incident_edges(network, connected_nodes3, mode = "in") 
    connected_edges2 <- incident_edges(network, connected_nodes5, mode = "out") 
    E(network)[unlist(connected_edges)]$color <- highlight_edge_color
    E(network)[unlist(connected_edges2)]$color <- highlight_edge_color
  } else {
    connected_edges <- incident_edges(network, all_upstream_nodes, mode = "in") 
    E(network)[unlist(connected_edges)]$color <- highlight_edge_color
  }
  
  
  # 设置相关节点和边的颜色
  #E(network)[unlist(connected_edges)]$color <- highlight_edge_color
  
  #设置网络图的布局
  if (is.null(netlayout)) {
    netlayout <- layout_(network, 
                         with_kk())
  }
  
  # 设置 PDF 的输出文件名、高度和宽度
  pdf.name <- paste0("ceRNA_highlight_", highlight_node, ".pdf")
  
  pdf(file = pdf.name, height = pdf_height, width = pdf_width)
  
  # 绘制网络图
  plot(network, edge.arrow.size = 0.4,
       layout = netlayout, 
       vertex.label.family = vertex_label_family, 
       res = c(1200, 1200))
  
  # 关闭 PDF 设备
  dev.off()
}
