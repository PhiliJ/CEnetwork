% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot_network.R
\name{plot_CEnetwork}
\alias{plot_CEnetwork}
\title{Plot ceRNA network}
\usage{
plot_CEnetwork(
  net,
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
  vertex.label.family = "sans"
)
}
\arguments{
\item{net}{A data frame that contains information about the nodes in the network. It should have two columns: "From" and "To". Use "UPnetwork", "DOWNnetwork" or "mynetwork".}

\item{nodeshape}{The shape of the nodes. Default is "circle". If you want to draw a 3D plot, use "sphere" instead.}

\item{circ_color}{a character specifying the color for circular RNA nodes}

\item{lnc_color}{a character specifying the color for long non-coding RNA nodes}

\item{mir_color}{a character specifying the color for microRNA nodes}

\item{Up_mrna_color}{a character specifying the color for upregulated mRNA nodes}

\item{Down_mrna_color}{a character specifying the color for downregulated mRNA nodes}

\item{circ_alpha}{a numeric value specifying the alpha transparency for circular RNA nodes}

\item{lnc_alpha}{a numeric value specifying the alpha transparency for long non-coding RNA nodes}

\item{mir_alpha}{a numeric value specifying the alpha transparency for microRNA nodes}

\item{Up_alpha}{a numeric value specifying the alpha transparency for upregulated mRNA nodes}

\item{Down_alpha}{a numeric value specifying the alpha transparency for downregulated mRNA nodes}

\item{circ_size}{a numeric value specifying the size for circular RNA nodes}

\item{lnc_size}{a numeric value specifying the size for long non-coding RNA nodes}

\item{mir_size}{a numeric value specifying the size for microRNA nodes}

\item{mrna_size}{a numeric value specifying the size for mRNA nodes}

\item{circ_labelsize}{a numeric value specifying the size for labels of circular RNA nodes}

\item{lnc_labelsize}{a numeric value specifying the size for labels of long non-coding RNA nodes}

\item{mir_labelsize}{a numeric value specifying the size for labels of microRNA nodes}

\item{mrna_labelsize}{a numeric value specifying the size for labels of mRNA nodes}

\item{label_color}{The color of the node labels. Default is "black".}

\item{edge_color}{The color of the edges. Default is "grey".}

\item{edge_alpha}{The alpha value of the edges. Default is 0.75.}

\item{netlayout}{The layout of the network.  Default is NULL, and the function will calculate the layout using the Kamada-Kawai algorithm.}

\item{pdf.name}{The name of the PDF file to be saved. Default is "TFnet.pdf".}

\item{pdf.height}{The height of the PDF file. Default is 15.}

\item{pdf.width}{The width of the PDF file. Default is 15.}

\item{edge.arrow.size}{The size of the arrows on the edges. Default is 0.4.}

\item{vertex.label.family}{The font family of the node labels. Default is "sans".}
}
\value{
A PDF file containing the network plot, and the network object stored in the global environment
}
\description{
This function plots a ceRNA network, including circRNAs, lncRNAs, miRNAs and mRNAs
}
