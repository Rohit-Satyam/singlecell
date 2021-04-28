## Steps for QCing single cell data

```r
suppressMessages(library("DropletUtils"))
suppressMessages(library(scater))
suppressMessages(library(SummarizedExperiment))
library(ggplot2)
library(scales)
```
## Setting up the data
We load in the raw count matrix using the `read10xCounts()` function from the _[DropletUtils](https://bioconductor.org/packages/3.10/DropletUtils)_ package. This will create a `SingleCellExperiment` object where each column corresponds to a cell barcode.
```r
setwd("~/rohit/kaust_proj")
root <- "~/rohit/kaust_proj"
p16D <- paste0(root,"/16D/raw_feature_bc_matrix/")  
p16R <- paste0(root,"/16R/raw_feature_bc_matrix/") 
p40D <- paste0(root,"/40D/raw_feature_bc_matrix/")  
p40R <- paste0(root,"/40R/raw_feature_bc_matrix/") 

p16D_sce <- read10xCounts(p16D, col.names=TRUE) #dim: 39250 6794880
p16R_sce <- read10xCounts(p16R, col.names=TRUE) #dim: 39250 6794880
p40D_sce <- read10xCounts(p40D, col.names=TRUE) #dim: 39250 6794880
p40R_sce <- read10xCounts(p40R, col.names=TRUE) #dim: 39250 6794880

#head(rowData(p16D_sce))
```
## Filtering for genes (Rows)
> Since the data contains both human and plasmodium genes, we will remove human genes because RBC got no organelle or DNA. Also we will remove genes that are not expressed in any cell: 
```r
name <- rownames(p16D_sce)
hg <- grep("hg38", rownames(p16D_sce))
name <- name[-hg]

## keep all the plasmodium genes
sce <-sce[name,]

p16D_sce <- p16D_sce[name,] #dim: 39250 6794880
p16R_sce <- p16R_sce[name,]
p40D_sce <- p40D_sce[name,]
p40R_sce <- p40R_sce[name,]


keep_feature <- rowSums(counts(p16D_sce) > 0) > 0
p16D_sce <- p16D_sce[keep_feature, ] #dim: 5446 6794880 

keep_feature <- rowSums(counts(p16R_sce) > 0) > 0
p16R_sce <- p16R_sce[keep_feature, ] #dim: 5408 6794880

keep_feature <- rowSums(counts(p40D_sce) > 0) > 0
p40D_sce <- p40D_sce[keep_feature, ] #dim: 5405 6794880

keep_feature <- rowSums(counts(p40R_sce) > 0) > 0
p40R_sce <- p40R_sce[keep_feature, ] #dim: 5363 6794880

## Re-lable the row names for readability. In case of Plasmodium the ID and Symbol column have same entry so it can be skipped: rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol). Use for human data.
```
## Calling for empty cells (Columns)
This is not entirely straightforward as empty droplets can contain ambient (i.e., extracellular) RNA that can be captured and sequenced. 
The waterfall plot shows the log-count against the log-rank of each barcode.The barcodes are ranked based on the number of count each barcode has.
We will compute these statistics using the `barcodeRanks` function from the `DropletUtils` package. And then we will extra the statistics and perform some plotting using the ggplot2 package.

Here given is a waterfall plot function that ingests the bcrank object from `bcrank <- barcodeRanks(counts(p40R_sce))`

Usage:
`plot_waterfall(bcrank,"p40R_sce")`

```r  
plot_waterfall <- function (bcrank, sample){ 
barcode_data = as.data.frame(bcrank)
barcode_points = data.frame(
  type = c("inflection", "knee"),
  value = c(bcrank@metadata$inflection, bcrank@metadata$knee))


ggplot(data = barcode_data, aes(x = rank, y = total)) +
  geom_point() +
  geom_hline(data = barcode_points,
             aes(yintercept = value,
                 colour = type), linetype = 2) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(title = paste0("Waterfall plot of read counts (log) for ", sample),
       x = "Log Rank",
       y = "Log counts")+theme_bw()+theme(axis.text.x = element_text(face="bold", color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=12))
       }
```
Making Waterfall plot for all 4 samples:

The waterfall plots show two points of interest:

- **inflection** point is the point on the curve where the first derivative is minimised
- **knee point** is the point where the signed curvature is minimised <br>
These two points potentially indicate `empty droplets` (background barcodes). To elaborate on this concept, it is helpful to know that for microfluidic sequencing systems each droplet contains a bead and each bead are marked uniquely with a barcode. Intuitively, we can consider a droplet and the barcode associated with it as being a unique cell. Due to the technology, however, sometimes the droplet don’t always contain a cell, or it might contain multiples cells. We can now use the waterfall plot to detect those empty droplets. The barcodes with the number of counts dropping rapidly beyond knee point may suggest empty reads.

>High quality barcodes are located on the left hand side of the plot, and thresholding is performed by identifying the “knee” on the curve. On the right hand side, past the inflection point, are barcodes which have relatively low numbers of reads, and are therefore considered to have had failure in capture and to be too noisy for further analysis. [Source](https://liorpachter.wordpress.com/tag/knee-plot/#:~:text=A%20single%2Dcell%20RNA%2Dseq,%E2%80%9Cknee%E2%80%9D%20on%20the%20curve.)

```r
bcrank <- barcodeRanks(counts(p16R_sce))
plot_waterfall(bcrank,"p16R_sce")

bcrank <- barcodeRanks(counts(p16D_sce))
plot_waterfall(bcrank,"p16D_sce")

bcrank <- barcodeRanks(counts(p40R_sce))
plot_waterfall(bcrank,"p40R_sce")

bcrank <- barcodeRanks(counts(p40D_sce))
plot_waterfall(bcrank,"p40D_sce")
```


## This gives us total UMI count for each barcode in the PBMC dataset, plotted against its rank (in decreasing order of total counts)

## emptyDrops() function to test whether the expression profile for each cell barcode is significantly different from the ambient RNA pool.

set.seed(100)
e.out <- emptyDrops(counts(sce))
is.cell <- sum(e.out$FDR <= 0.001, na.rm=TRUE)

## dim before filtering:dim: 5712 6794880 and after filtering dim: 5712 4255. Cross-checked with the filterred counts too iif same number of cells are rescued.

plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")
##Droplets detected as cells should show up with large negative log-probabilities or very large total counts (based on the knee point reported by barcodeRanks
```

## Quality control on the cells


```{r, eval=FALSE}
is.mito <- grep("mito", rownames(sce))
#per.cell <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
sce <- calculateQCMetrics(sce, feature_controls=list(Mito=is.mito))
#summary(per.cell$sum)
par(mfrow=c(1,3))
hist(sce$log10_total_counts, breaks=20, col="grey80",
    xlab="Log-total UMI count")
hist(sce$log10_total_features_by_counts, breaks=20, col="grey80",
    xlab="Log-total number of expressed features")
hist(sce$pct_counts_Mito, breaks=20, col="grey80",
    xlab="Proportion of reads in mitochondrial genes")

#using less stringent criteria to filter cells with higher mitochondrial genes.
high.mito <- isOutlier(sce$pct_counts_Mito, nmads=3, type="higher")
sce <- sce[,!high.mito]
summary(high.mito)

ave <- calculateAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80")

#most highly expressed genes is dominated
plotHighestExprs(sce, n=25) + theme(text = element_text(size=14))
```
