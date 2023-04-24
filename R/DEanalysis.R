#' @title Construct ceRNA network that targets your mRNA of interest
#' @description Construct ceRNA network that targets your mRNA of interest (only support 1 mRNA, I will update later)
#' @param file a character string specifying the name of the input file to be read
#' @param min_expr a numeric value specifying the minimum expression level to consider a gene expressed (default: 1)
#' @param min_samples a numeric value specifying the minimum number of samples where a gene should be expressed (default: 1)
#' @importFrom base apply rowSums
#' @return a filtered gene expression data frame with genes as rows and samples as columns, where only genes with expression levels greater than or equal to min_expr in at least min_samples samples are retained. The filtered data frame is also written to a file named "<input_file_name>_filtered.txt".
#' @export
readexp <- function(file, min_expr = 1, min_samples = 1){
  data <- read.table(file, header = TRUE, row.names = NULL)
  data <- data[!duplicated(data$SYMBOL),]
  rownames(data) <- data$SYMBOL
  data <- data[, -1]
  data <- data[rowSums(data)>min_samples,]
  data <- data[apply(data,1, function(x) sum(x>min_expr) > min_samples),] 
  new_file <- paste0(sub(".txt", "", file), "_filtered.txt")
  write.table(data, file = new_file, sep = "\t")
  return(data)
}


#' @title DEG Analysis of RNA-seq Data using DESeq with P value
#' @description This function performs differential expression analysis using the DESeq2 package.
#' @param counts A count matrix where rows are genes and columns are samples.
#' @param phenodata A data frame containing the phenotype data for each sample.
#' @param complist A character vector containing the names of the samples to compare.
#' @param logfc The log2 fold change cutoff for differential expression. Default is 0.
#' @param pval The adjusted p-value or p-value cutoff for differential expression. Default is 0.05.
#' @param padj Logical indicating whether to use adjusted p-values or not. Default is TRUE.
#' @import DESeq2
#' @return DEresult: A list containing the following items:
#' table: A data frame with the results of the differential expression analysis.
#' rank_vector: A vector containing the rank of each gene based on the effect size.
#' degs: A character vector containing the names of the differentially expressed genes.
#' e: A matrix of normalized expression values.
#' up: A character vector containing the names of the upregulated genes.
#' down: A character vector containing the names of the downregulated genes.
#' nosig: A character vector containing the names of the non-differentially expressed genes.
#' @export
#' @import DESeq2
DiffDESeq2 = function(counts, phenodata, complist,logfc = 0, pval = 0.05, padj = TRUE ) {
  
  phenodata = data.frame(celltype = factor(group), row.names = samplenames)
  dds = DESeqDataSetFromMatrix(countData = counts, colData = phenodata, design = ~ celltype)
  dds <- estimateSizeFactors(dds)
  e <- counts(dds, normalized = TRUE)
  
  dds = DESeq(dds)
  if (padj) {
    DEtable = results(dds, alpha = pval)
    DEtable = data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$padj), ]
    degs = rownames(DEtable[abs(DEtable$log2FoldChange) > logfc & DEtable$padj < pval, ])
  } else {
    DEtable = results(dds, alpha = pval)
    DEtable = data.frame(DEtable[complete.cases(DEtable), ])
    DEtable = DEtable[order(DEtable$pvalue), ]
    degs = rownames(DEtable[abs(DEtable$log2FoldChange) > logfc & DEtable$pvalue < pval, ])
  }
  
  # 计算高低差异基因数目
  upregulated = rownames(DEtable[DEtable$log2FoldChange > logfc & (if (padj) DEtable$padj < pval else DEtable$pvalue < pval), ])
  downregulated = rownames(DEtable[DEtable$log2FoldChange < -logfc & (if (padj) DEtable$padj < pval else DEtable$pvalue < pval), ])
  no_difference = rownames(DEtable[abs(DEtable$log2FoldChange) <= logfc | (if (padj) DEtable$padj >= pval else DEtable$pvalue >= pval), ])
  
  # 输出结果
  num_upregulated = length(upregulated)
  num_downregulated = length(downregulated)
  num_no_difference = length(no_difference)
  
  print(paste("Down:", num_downregulated))
  print(paste("NotSig:", num_no_difference))
  print(paste("Up:", num_upregulated))
  
  rank_vector = abs(DEtable$stat); names(rank_vector) = rownames(DEtable)
  
  return(list(table = DEtable, rank_vector = rank_vector, degs = degs, e = e,
              up = upregulated, down = downregulated, nosig = no_difference))
}










