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

# split the dataset into a list of two seurat objects (stim and CTRL)
sample.list <- SplitObject(seurat, split.by = "Sample")

sample.list <- lapply(X = sample.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = sample.list)
plasmodium.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = features)
plasmodium.combined <- IntegrateData(anchorset = plasmodium.anchors)
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(plasmodium.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
plasmodium.combined <- ScaleData(plasmodium.combined, verbose = FALSE)
plasmodium.combined <- RunPCA(plasmodium.combined, npcs = 30, verbose = FALSE)

### Running Harmony

pf_hmony <- plasmodium.combined %>% RunHarmony("Sample", plot_convergence = TRUE)








plasmodium.combined <- RunUMAP(plasmodium.combined, reduction = "pca", dims = 1:30)
plasmodium.combined <- FindNeighbors(plasmodium.combined, reduction = "pca", dims = 1:30)
plasmodium.combined <- FindClusters(plasmodium.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(plasmodium.combined, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(plasmodium.combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(plasmodium.combined, reduction = "umap", group.by = "batch")
p1 + p2 + p3

```
