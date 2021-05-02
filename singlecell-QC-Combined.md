

```r
library(DropletUtils)
seed <- 100
out_path <- here::here()
data_path <- "/home/pichkari/rohit/kaust_proj/"
## I moved the gz files from raw_data folder to samplefolder
dirs <- c("16D", "16R", "40D", "40R")
sce <- DropletUtils::read10xCounts(samples=paste0(data_path,dirs), col.names = TRUE)
```
# Add metadata
```r
meta <- read_excel(paste0(data_path, "EXP_CR005 Samples ID.xlsx"))
sce$Sample <- gsub(".*kaust_proj/","", sce$Sample)
t <- plyr::count(sce$Sample)
count <- unique(t$freq)
sce$batch <- c(rep("16hpi",count),rep("16hpi",count),rep("40hpi",count),rep("40hpi",count))
colnames(sce) <- paste0(sce$batch, ".", sce$Barcode)
#rownames(sce) <- paste0(rowData(sce)$ID, ".", rowData(sce)$Symbol)

sce$batch <- factor(sce$batch) 
dim(sce) #dim: 39250 27179520
```

# Filtering
```r
keep_features <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_features, ] #dim: 6900 27179520

## Removing human genes
name <- rownames(sce)
hg <- grep("hg38", rownames(sce))
name <- name[-hg]
sce <-sce[name,]#dim: 5543 27179520
```
# Diagnosing for empty cells (Columns)

```r
## Making quick QC plots

bcrank <- barcodeRanks(counts(sce))
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
       
       
       plot_diag <- function(sceobj){
e.out <- emptyDrops(counts(sceobj))
is.cell <- sum(e.out$FDR <= 0.001, na.rm=TRUE)
plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"),
     xlab="Total UMI count", ylab="-Log Probability")
     }

empty_filt <- function(sceobj){
e.out <- emptyDrops(counts(sceobj))
sceobj <- sceobj[,which(e.out$FDR <= 0.001)]
return(sceobj)}

sce <- empty_filt(sce)
dim(sce)

#[1]  5543 19227

table(sce$batch)

#16hpi 40hpi 
#10740  8487 

table(sce$Sample)

# 16D  16R  40D  40R 
#4779 5961 4361 4126

```
# Percell QC

```r
is.mito <- grep("mito", rownames(sce))
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito))

hist(sce$subsets_Mito_percent, breaks=20, col="grey80",
     xlab="Proportion of reads in mitochondrial genes")
hist(sce$total, breaks=20, col="grey80",
     xlab="Log-total number of counts for the cell (i.e., the library size)")
     ave <- calculateAverage(sce)
rowData(sce)$AveCount <- ave
hist(log10(ave), col="grey80")
```
# Normalising for cell Biases

```r
library(scran)
library(BiocSingular)
set.seed(1000)
clusters <- quickCluster(sce, BSPARAM=IrlbaParam())
table(clusters)

sce <- computeSumFactors(sce, min.mean=0.1, cluster=clusters)
summary(sizeFactors(sce))

sce <- logNormCounts(sce)
dec <- modelGeneVarByPoisson(sce)
plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
top.hvgs <- getTopHVGs(dec, prop=0.1)
# Evaluate PCs
sce2 <- denoisePCA(sce, subset.row = top.hvgs, technical = dec, BSPARAM=IrlbaParam())
# make TSNE plot
sce2 <- runTSNE(sce2, dimred = "PCA")
# make UMAP plot
sce2 <- runUMAP(sce2, dimred = "PCA")

g <- buildSNNGraph(sce2, k = 10, use.dimred = "PCA")
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce2) <- factor(clust)
# 
colData(sce2)
plotUMAP(sce2, colour_by = "label")
plotTSNE(sce2, colour_by = "label")
top.dec <- subset(dec, rownames(dec) %in% top.hvgs)
plotExpression(sce2, features=rownames(top.dec)[1:10])

## Number of clusters can be seen by
table(sce2$label)
```

# Remove doublets

```r
dbl.dens <- scDblFinder::computeDoubletDensity(sce2,  d=ncol(reducedDim(sce2)),subset.row=top.hvgs)
sce2$DoubletScore <- dbl.dens
plotUMAP(sce2, colour_by = "DoubletScore")
## In this analysis doublets are rare.

#We can clean our data up based on this 95% quantile cut-off.
plotColData(sce2, x = "label", y = "DoubletScore", colour_by = "label") + geom_hline(yintercept = quantile(colData(sce2)$DoubletScore, 0.95), lty = "dashed", color = "red")
cut_off <- quantile(sce2$DoubletScore, 0.95)
sce2$isDoublet <- c("no", "yes")[factor(as.integer(sce2$DoubletScore >= cut_off),levels = c(0, 1))]
table(sce2$isDoublet)

#   no   yes 
#18259   968
```
