## Finding the Var genes:

We downloaded the gff file from PlasmoDB v37

```bash
grep PfEMP1 PlasmoDB-37_Pfalciparum3D7.gff | awk '{print $9}' | sort -u | cut -f '2' -d'=' | cut -f '1' -d';' | cut -f '1' -d'.' | sort -u | wc -l > var_genes.csv
```

This gave us 105 var genes

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

##pf_hmony <- plasmodium.combined %>% RunHarmony("Sample", plot_convergence = TRUE)








plasmodium.combined <- RunUMAP(plasmodium.combined, reduction = "pca", dims = 1:30)
plasmodium.combined <- FindNeighbors(plasmodium.combined, reduction = "pca", dims = 1:30)
plasmodium.combined <- FindClusters(plasmodium.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(plasmodium.combined, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(plasmodium.combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(plasmodium.combined, reduction = "umap", group.by = "batch")
p1 + p2 + p3

```

## Integrating with the malaria dataset

```r
## I am told that this dataset of Malaria Cell Atlas is already preprocessed so we need not to preprocess it again. Also it doesnot contain mitochondrial genes
dat.m = read.csv("~/Desktop/Plasmodium_scRNAseq/malaria cell atlas data/pf10xIDC_counts.csv", header = TRUE ,row.names = 1)
sobj.m <- CreateSeuratObject(counts = dat.m, project = "Malaria-Cell-Atlas", min.cells = 5, min.features = 100)
sobj.m$condition = "malaria.atlas"
sobj.m[["percent.mt"]] <- PercentageFeatureSet(object = sobj.m, pattern = "mito")

## Loading the sce_clean data obtained from preprocessed steps 
counts <- assays(sce_clean)[[1]]

# Remember to check the uniformity of row features while doing such integrative analysis. Here we remove these extra characters that we had in gene names. Since you can't change the gene names post creation of seurat object, better to do it before creating one: Source: https://github.com/satijalab/seurat/issues/2617
rownames(counts) <- gsub("pf3d7_","",rownames(counts))
sobj.ie <- CreateSeuratObject(counts = counts, project = "Inlab-Experiment", min.cells = 5, min.features = 100)
t2 <- sobj.ie@meta.data
t <- data.frame(colData(sce_clean))
t3 <- t[rownames(t2),]
sobj.ie@meta.data <- t3

## 
set.seed(1000)
sample.list <- SplitObject(sobj.ie, split.by = "Sample")
sample.list$mca = sobj.m
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
plasmodium.combined <- RunUMAP(plasmodium.combined, reduction = "pca", dims = 1:30)
plasmodium.combined <- FindNeighbors(plasmodium.combined, reduction = "pca", dims = 1:30)
plasmodium.combined <- FindClusters(plasmodium.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(plasmodium.combined, reduction = "umap", group.by = "Sample")
p2 <- DimPlot(plasmodium.combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(plasmodium.combined, reduction = "umap", group.by = "batch")
p1+p2+p3
```

