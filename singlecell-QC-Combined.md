

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

```
       

