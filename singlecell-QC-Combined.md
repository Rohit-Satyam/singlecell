

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



