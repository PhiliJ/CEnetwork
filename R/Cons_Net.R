#' @title Construct ceRNA network that targets up-regulated mRNAs
#' @description Construct ceRNA network that targets up-regulated mRNAs
#' @param clipExpNum_mir2m_cutoff the minimum CLIP-seq experiment number for miRNA-mRNA interactions (default = 1)
#' @param databaseNum_mir2m_cutoff the minimum number of databases supporting miRNA-mRNA interactions (default = 3)
#' @param clipExpNum_lnc2mir_cutoff the minimum CLIP-seq experiment number for lncRNA-miRNA interactions (default = 1)
#' @param clipExpNum_circ2mir_cutoff the minimum CLIP-seq experiment number for circRNA-miRNA interactions (default = 1)
#'
#' @import dplyr
#' @import tidyverse
#' @import tidyr
#' @importFrom datasets miRNA_mRNA_mouse lncRNA_miRNA_mouse circRNA_miRNA_mouse
#' @return UPnetwork: ceRNA network that targets up-regulated mRNAs
#' @export 
Cons_Up_net <- function(clipExpNum_mir2m_cutoff = 1,
                        databaseNum_mir2m_cutoff = 3,
                        clipExpNum_lnc2mir_cutoff = 1, 
                        clipExpNum_circ2mir_cutoff = 1
){
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  load("miRNA_mRNA_mouse.rdata")
  load("lncRNA_miRNA_mouse.rdata")
  load("circRNA_miRNA_mouse.rdata")
  
  sub <- subset(miRNA_mRNA_mouse, clipExpNum >= clipExpNum_mir2m_cutoff & (PITA + RNA22 + miRmap + microT + miRanda + PicTar + TargetScan) >= databaseNum_mir2m_cutoff)
  
  UPm2mir <- sub[which(sub$mRNA %in% DEm$up),]
  UPm2mir <-  UPm2mir[, c("miRNA", "mRNA")]
  UPm2Dmir <- UPm2mir[which(UPm2mir$miRNA %in% DEmir$down),]
  
  sub2 <- subset(lncRNA_miRNA_mouse, clipExpNum >= clipExpNum_lnc2mir_cutoff)
  
  Dmir2lnc <- sub2[which(sub2$miRNA %in% UPm2Dmir$miRNA),]
  Dmir2Hlnc <- Dmir2lnc[which(Dmir2lnc$lncRNA %in% DElnc$up),]
  
  sub3 <- subset(circRNA_miRNA_mouse, clipExpNum >= clipExpNum_circ2mir_cutoff)
  
  Dmir2circ <- sub3[which(sub3$miRNA %in% UPm2Dmir$miRNA),]
  Dmir2Hcirc <- Dmir2circ[which(Dmir2circ$circRNA %in% DEcirc$up),]
  
  UPm2Dmir <- subset(UPm2mir, miRNA %in% c(Dmir2Hlnc$miRNA, Dmir2Hcirc$miRNA))
  
  #Construct Network of UP mRNAs
  mir2m <- UPm2Dmir[, c("miRNA", "mRNA")]
  colnames(mir2m) <- c("From", "To")
  lnc2mir <- Dmir2Hlnc[, c("lncRNA", "miRNA")]
  colnames(lnc2mir) <- c("From", "To")
  circ2mir <- Dmir2Hcirc[, c("circRNA", "miRNA")]
  colnames(circ2mir) <- c("From", "To")
  
  UPnetwork <- rbind(mir2m, lnc2mir, circ2mir)
  
  write.table(UPnetwork, file = "UPnetwork.txt", sep = "\t", row.names = FALSE)
  return(UPnetwork)
  
}


#' @title Construct ceRNA network that targets down-regulated mRNAs
#' @description Construct ceRNA network that targets down-regulated mRNAs
#' @param clipExpNum_mir2m_cutoff the minimum CLIP-seq experiment number for miRNA-mRNA interactions (default = 1)
#' @param databaseNum_mir2m_cutoff the minimum number of databases supporting miRNA-mRNA interactions (default = 3)
#' @param clipExpNum_lnc2mir_cutoff the minimum CLIP-seq experiment number for lncRNA-miRNA interactions (default = 1)
#' @param clipExpNum_circ2mir_cutoff the minimum CLIP-seq experiment number for circRNA-miRNA interactions (default = 1)
#'
#' @import dplyr
#' @import tidyverse
#' @import tidyr
#' @importFrom datasets miRNA_mRNA_mouse lncRNA_miRNA_mouse circRNA_miRNA_mouse
#' @return DOWNnetwork: ceRNA network that targets down-regulated mRNAs
#' @export 
Cons_Down_net <- function(clipExpNum_mir2m_cutoff = 1,
                          databaseNum_mir2m_cutoff = 3,
                          clipExpNum_lnc2mir_cutoff = 1, 
                          clipExpNum_circ2mir_cutoff = 1
){
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  load("miRNA_mRNA_mouse.rdata")
  load("lncRNA_miRNA_mouse.rdata")
  load("circRNA_miRNA_mouse.rdata")
  
  sub <- subset(miRNA_mRNA_mouse, clipExpNum >= clipExpNum_mir2m_cutoff & (PITA + RNA22 + miRmap + microT + miRanda + PicTar + TargetScan) >= databaseNum_mir2m_cutoff)
  
  Dm2mir <- sub[which(sub$mRNA %in% DEm$down),]
  Dm2mir <-  Dm2mir[, c("miRNA", "mRNA")]
  Dm2Hmir <- Dm2mir[which(Dm2mir$miRNA %in% DEmir$up),]
  
  sub2 <- subset(lncRNA_miRNA_mouse, clipExpNum >= clipExpNum_lnc2mir_cutoff)
  
  Hmir2lnc <- sub2[which(sub2$miRNA %in% Dm2Hmir$miRNA),]
  Hmir2Dlnc <- Hmir2lnc[which(Hmir2lnc$lncRNA %in% DElnc$down),]
  
  sub3 <- subset(circRNA_miRNA_mouse, clipExpNum >= clipExpNum_circ2mir_cutoff)
  
  Hmir2circ <- sub3[which(sub3$miRNA %in% Dm2Hmir$miRNA),]
  Hmir2Dcirc <- Hmir2circ[which(Hmir2circ$circRNA %in% DEcirc$down),]
  
  Dm2Hmir <- subset(Dm2Hmir, miRNA %in% c(Hmir2Dlnc$miRNA, Hmir2Dcirc$miRNA))
  
  #Construct Network of UP mRNAs
  mir2m <- Dm2Hmir[, c("miRNA", "mRNA")]
  colnames(mir2m) <- c("From", "To")
  lnc2mir <- Hmir2Dlnc[, c("lncRNA", "miRNA")]
  colnames(lnc2mir) <- c("From", "To")
  circ2mir <- Hmir2Dcirc[, c("circRNA", "miRNA")]
  colnames(circ2mir) <- c("From", "To")
  
  DOWNnetwork <- rbind(mir2m, lnc2mir, circ2mir)
  
  write.table(DOWNnetwork, file = "DOWNnetwork.txt", sep = "\t", row.names = FALSE)
  return(DOWNnetwork)
  
}

#' @title Construct ceRNA network that targets your mRNA of interest
#' @description Construct ceRNA network that targets your mRNA of interest (only support 1 mRNA, I will update later)
#' @param my_mRNA Your mRNA of interest (only support 1 mRNA, I will update later)
#' @param clipExpNum_mir2m_cutoff the minimum CLIP-seq experiment number for miRNA-mRNA interactions (default = 1)
#' @param databaseNum_mir2m_cutoff the minimum number of databases supporting miRNA-mRNA interactions (default = 3)
#' @param clipExpNum_lnc2mir_cutoff the minimum CLIP-seq experiment number for lncRNA-miRNA interactions (default = 1)
#' @param clipExpNum_circ2mir_cutoff the minimum CLIP-seq experiment number for circRNA-miRNA interactions (default = 1)
#'
#' @import dplyr
#' @import tidyverse
#' @import tidyr
#' @importFrom datasets miRNA_mRNA_mouse lncRNA_miRNA_mouse circRNA_miRNA_mouse
#' @return mynetwork: ceRNA network that targets your mRNA of interest
#' @export 
Cons_my_net <-  function(my_mRNA,
                         clipExpNum_mir2m_cutoff = 1,
                         databaseNum_mir2m_cutoff = 3,
                         clipExpNum_lnc2mir_cutoff = 1, 
                         clipExpNum_circ2mir_cutoff = 1)
{ 
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  load("miRNA_mRNA_mouse.rdata")
  load("lncRNA_miRNA_mouse.rdata")
  load("circRNA_miRNA_mouse.rdata")
  
  if (my_mRNA %in% DEm$up) {
    sub <- subset(miRNA_mRNA_mouse, clipExpNum >= clipExpNum_mir2m_cutoff & (PITA + RNA22 + miRmap + microT + miRanda + PicTar + TargetScan) >= databaseNum_mir2m_cutoff)
    
    UPm2mir <- sub[which(sub$mRNA %in% my_mRNA),]
    UPm2mir <-  UPm2mir[, c("miRNA", "mRNA")]
    UPm2Dmir <- UPm2mir[which(UPm2mir$miRNA %in% DEmir$down),]
    
    sub2 <- subset(lncRNA_miRNA_mouse, clipExpNum >= clipExpNum_lnc2mir_cutoff)
    
    Dmir2lnc <- sub2[which(sub2$miRNA %in% UPm2Dmir$miRNA),]
    Dmir2Hlnc <- Dmir2lnc[which(Dmir2lnc$lncRNA %in% DElnc$up),]
    
    sub3 <- subset(circRNA_miRNA_mouse, clipExpNum >= clipExpNum_circ2mir_cutoff)
    
    Dmir2circ <- sub3[which(sub3$miRNA %in% UPm2Dmir$miRNA),]
    Dmir2Hcirc <- Dmir2circ[which(Dmir2circ$circRNA %in% DEcirc$up),]
    
    UPm2Dmir <- subset(UPm2mir, miRNA %in% c(Dmir2Hlnc$miRNA, Dmir2Hcirc$miRNA))
    
    #Construct Network of UP mRNAs
    mir2m <- UPm2Dmir[, c("miRNA", "mRNA")]
    colnames(mir2m) <- c("From", "To")
    lnc2mir <- Dmir2Hlnc[, c("lncRNA", "miRNA")]
    colnames(lnc2mir) <- c("From", "To")
    circ2mir <- Dmir2Hcirc[, c("circRNA", "miRNA")]
    colnames(circ2mir) <- c("From", "To")
    
    mynetwork <- rbind(mir2m, lnc2mir, circ2mir)
    
    write.table(mynetwork, file = "mynetwork.txt", sep = "\t", row.names = FALSE)
    return(mynetwork)
    
  } else {
    sub <- subset(miRNA_mRNA_mouse, clipExpNum >= clipExpNum_mir2m_cutoff & (PITA + RNA22 + miRmap + microT + miRanda + PicTar + TargetScan) >= databaseNum_mir2m_cutoff)
    
    Dm2mir <- sub[which(sub$mRNA %in% my_mRNA),]
    Dm2mir <-  Dm2mir[, c("miRNA", "mRNA")]
    Dm2Hmir <- Dm2mir[which(Dm2mir$miRNA %in% DEmir$up),]
    
    sub2 <- subset(lncRNA_miRNA_mouse, clipExpNum >= clipExpNum_lnc2mir_cutoff)
    
    Hmir2lnc <- sub2[which(sub2$miRNA %in% Dm2Hmir$miRNA),]
    Hmir2Dlnc <- Hmir2lnc[which(Hmir2lnc$lncRNA %in% DElnc$down),]
    
    sub3 <- subset(circRNA_miRNA_mouse, clipExpNum >= clipExpNum_circ2mir_cutoff)
    
    Hmir2circ <- sub3[which(sub3$miRNA %in% Dm2Hmir$miRNA),]
    Hmir2Dcirc <- Hmir2circ[which(Hmir2circ$circRNA %in% DEcirc$down),]
    
    Dm2Hmir <- subset(Dm2Hmir, miRNA %in% c(Hmir2Dlnc$miRNA, Hmir2Dcirc$miRNA))
    
    #Construct Network of UP mRNAs
    mir2m <- Dm2Hmir[, c("miRNA", "mRNA")]
    colnames(mir2m) <- c("From", "To")
    lnc2mir <- Hmir2Dlnc[, c("lncRNA", "miRNA")]
    colnames(lnc2mir) <- c("From", "To")
    circ2mir <- Hmir2Dcirc[, c("circRNA", "miRNA")]
    colnames(circ2mir) <- c("From", "To")
    
    mynetwork <- rbind(mir2m, lnc2mir, circ2mir)
    
    write.table(mynetwork, file = paste0("mynetwork_", my_mRNA, ".txt"), sep = "\t", row.names = FALSE)
    return(mynetwork)
    
  }
} 


