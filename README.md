# CEnetwork
CEnetwork: an R package that constructs and plots ceRNA network based on RNA-seq data and ENCORI database.

`CEnetwork` is an R package for Differential-expression analysis (DE-analysis) using DESeq2 (will update limma in the future), construction of ceRNA network (based on DE-analysis and ENCORI database) and ceRNA network plotting, all in one. 

What you need to do is simple: just input the expression matrix of circRNA, lncRNA miRNA and mRNA, and you will get the differentially-expression genes(DEGs), ceRNA network that targets up-regulated mRNAs, down-regulated mRNAs or your mRNA of interest, a normal ceRNA network plot and network plot with highlighted RNAs of interest.

`CEnetwork` is licensed under the GNU General Public License v3.0 https://github.com/PhiliJ/TFnetwork/blob/main/LICENSE

### Author

Li Lei <2020320243@stu.cqmu.edu.cn>

### Installation:

You can install `CEnetwork` like so:

``` r
library(devtools)
install_github("PhiliJ/CEnetwork", dependencies = T)
``` 
### Step 1 Preparation of expreesion matrix for later analyses

For example, for treated group and NC group:

First thing first, you will need to prepare the expression table of circRNAs, lncRNAs, miRNAs and mRNAs. 
The first row of the table should be named as "SYMBOL", which are symbols of RNAs, and from the second row
should be the row count of each gene, then name these table like below. Maybe it's better that you DON'T change the name of the inputs.

``` r
library(CEnetwork)
setwd("your path to the expression files")

# read and preprocess expression matrix
lncRNA <- readexp("lncRNA_exp.txt")
mRNA <- readexp("mRNA_exp.txt")
circRNA <- readexp("circRNA_exp.txt")
miRNA <- readexp("miRNA_exp.txt")
```


### Step 2 DE-analysis for DEGs
You can use `DESeq2` for DE-analysis. I will update `limma` in the future.
#### 2.1 Set up expression matrix and phenodata
``` r
group <- c(rep("nc", 3), rep("treat", 3))
samplenames <- colnames(mRNA)[1:6]
compList <- c("treat-nc")
```

#### 2.2 DE-analysis using DESeq2
The output `_norm_exp.txt` is the normalized expression.
``` r
DElnc <- DiffDESeq2 (counts = lncRNA, complist = compList, 
                     logfc = 2, pval = 0.05, padj = TRUE)
DEm <- DiffDESeq2 (counts = mRNA, complist = compList, 
                   logfc = 10, pval = 0.00001, padj = TRUE)
DEmir <- DiffDESeq2 (counts = miRNA, complist = compList, 
                     logfc = 2, pval = 0.05, padj = TRUE)
DEcirc <- DiffDESeq2 (counts = circRNA, complist = compList, 
                      logfc = 1, pval = 0.05, padj = TRUE)

write.table(DElnc$table, file = "DESeq2_lncRNA_result.txt", sep = "\t")
write.table(DElnc$e, file = "DESeq2_lncRNA_norm_exp.txt", sep = "\t")
write.table(DEm$table, file = "DESeq2_mRNA_result.txt", sep = "\t")
write.table(DEm$e, file = "DESeq2_mRNA_norm_exp.txt", sep = "\t")
write.table(DEmir$table, file = "DESeq2_miRNA_result.txt", sep = "\t")
write.table(DEmir$e, file = "DESeq2_miRNA_norm_exp.txt", sep = "\t")
write.table(DEcirc$table, file = "DESeq2_circRNA_result.txt", sep = "\t")
write.table(DEcirc$e, file = "DESeq2_circRNA_norm_exp.txt", sep = "\t")
```
