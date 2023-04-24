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

### Step 3 Construct ceRNA network
For now we got 3 options. You can use `Cons_Up_net` to construct ceRNA network that targets up-regulated mRNAs, `Cons_Down_net`
to construct ceRNA network that targets down-regulated mRNAs, and `Cons_my_net` to construct ceRNA network that targets your mRNA
of interest. However, only one mRNA of interest is allowed for `Cons_my_net` since this is the beta version. I will update support
for multiple mRNAs of interest, even miRNA, lncRNA or circRNA in the near future.

You may use `my_mRNA =` to define your mRNA of interest. 

`Cons_Up_net`, `Cons_Down_net` and `Cons_my_net` will also return an '.txt' file contains source nodes and targets, 
and you can use this to generate network plot using other software like `Cytoscape`.
``` r
UPnetwork <- Cons_Up_net()
DOWNnetwork <- Cons_Down_net()
mynetwork <- Cons_my_net(my_mRNA = "Lancl3")
``` 



### Step 4 Plot ceRNA network
#### 4.1 Plot ceRNA network

Use `plot_CEnetwork`. You may need to check the details using `?plot_CEnetwork`.

Notebly, `plot_CEnetwork` add a constructed igraph object named "network" to the global environment of R, which is the basis of the next
section `plot_highlighted_CEnetwork`. So, you must draw the basic ceRNA network plot using  `plot_CEnetwork` in the first place, then 
`plot_highlighted_CEnetwork` will be functional.

Please Check `?plot_CEnetwork` for more details.

``` r
plot_CEnetwork(net = UPnetwork,
               pdf.name = "ceRNA_up.pdf")
#or
plot_CEnetwork(net = DOWNnetwork,
               pdf.name = "ceRNA_down.pdf")
#or
plot_CEnetwork(net = mynetwork,
               pdf.name = "ceRNA_my.pdf")
``` 

#### 4.2 Plot highlighted ceRNA network

Use `plot_highlighted_CEnetwork`. `CEnetwork` will automatically determine whether your node of interest is a upstream or downstream. You may need to check the details using `?plot_highlighted_CEnetwork`.

Notebly, `plot_highlighted_CEnetwork` is based on the igraph object named "network" (The network file constructed before using `plot_CEnetwork`). As a result, 
on one hand, the highlighted node must be in the "network", which means if you drawed the "UP" network, the node you want to highlight must be in the 'ceRNA_up.pdf', or `plot_highlighted_CEnetwork` will not be functional. On the other hand, like mentioned before, you must use `plot_CEnetwork` in the first place, then `plot_highlighted_CEnetwork` will be functional.


``` r
plot_highlighted_CEnetwork(network, highlight_node = "Snhg14")
plot_highlighted_CEnetwork(network, highlight_node = "mmu-miR-9-5p")
plot_highlighted_CEnetwork(network, highlight_node = "Cdh1")
```
[ceRNA_highlight_mmu-miR-9-5p.pdf](https://github.com/PhiliJ/CEnetwork/files/11305535/ceRNA_highlight_mmu-miR-9-5p.pdf)
