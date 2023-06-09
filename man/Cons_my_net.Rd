% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Cons_Net.R
\name{Cons_my_net}
\alias{Cons_my_net}
\title{Construct ceRNA network that targets your mRNA of interest}
\usage{
Cons_my_net(
  my_mRNA,
  clipExpNum_mir2m_cutoff = 1,
  databaseNum_mir2m_cutoff = 3,
  clipExpNum_lnc2mir_cutoff = 1,
  clipExpNum_circ2mir_cutoff = 1
)
}
\arguments{
\item{my_mRNA}{Your mRNA of interest (only support 1 mRNA, I will update later)}

\item{clipExpNum_mir2m_cutoff}{the minimum CLIP-seq experiment number for miRNA-mRNA interactions (default = 1)}

\item{databaseNum_mir2m_cutoff}{the minimum number of databases supporting miRNA-mRNA interactions (default = 3)}

\item{clipExpNum_lnc2mir_cutoff}{the minimum CLIP-seq experiment number for lncRNA-miRNA interactions (default = 1)}

\item{clipExpNum_circ2mir_cutoff}{the minimum CLIP-seq experiment number for circRNA-miRNA interactions (default = 1)}
}
\value{
mynetwork: ceRNA network that targets your mRNA of interest
}
\description{
Construct ceRNA network that targets your mRNA of interest (only support 1 mRNA, I will update later)
}
