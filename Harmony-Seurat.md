## Converting Single Cell Experiment object to Seurat Object

There is a direct method but that didn't work for me

```r
library(Seurat)

seurat <- as.Seurat(sce_clean)
```

### Indirect Way

```r
counts <- assays(sce_clean)[[1]]
seurat <- CreateSeuratObject(counts = counts, project = "Harmony_All", min.cells = 5)
t2 <- seurat@meta.data
t <- data.frame(colData(sce_clean))
t3 <- data.frame(t,t2)
seurat@meta.data <- t3

sample.list <- SplitObject(seurat, split.by = "batch")

sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

```
