---
title: "swapped_corrected"
author: "Rohit Satyam"
date: "06/05/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Data Preprocessing

```{r libraries}
suppressPackageStartupMessages({
 # library(scran) # giving error
  library(viridis)
  library(multtest)
  library(harmony)
  library(cowplot)
  library(scDblFinder)
library(BiocSingular)
    library("DropletUtils")
    library(scater)
    library(SummarizedExperiment)
    library(ggplot2)
    library(scales)
    library(Seurat)
  })

```

### Reading Data

```{r read}
set.seed(1000)
root <- here::here()  #"/home/pichkari/rohit/kaust_proj/sample_swap_corrected"
p16D <- paste0(root,"/16D/raw_feature_bc_matrix/")  
p16R <- paste0(root,"/16R/raw_feature_bc_matrix/") 
p40D <- paste0(root,"/40D/raw_feature_bc_matrix/")  
p40R <- paste0(root,"/40R/raw_feature_bc_matrix/") 
p16D_sce <- read10xCounts(p16D, col.names=TRUE) #dim: 39250 6794880
p16R_sce <- read10xCounts(p16R, col.names=TRUE) #dim: 39250 6794880
p40D_sce <- read10xCounts(p40D, col.names=TRUE) #dim: 39250 6794880
p40R_sce <- read10xCounts(p40R, col.names=TRUE) #dim: 39250 6794880
```

### Filtering for non-expressed genes

```{r}
p16D_sce <- filter(p16D_sce)
p16R_sce <- filter(p16R_sce)
p40D_sce <- filter(p40D_sce)
p40R_sce <- filter(p40R_sce)

```
### Diagnosing for empty cells (Columns)

```{r}
bcrank_p16D <- barcodeRanks(counts(p16D_sce))
bcrank_p16R <- barcodeRanks(counts(p16R_sce))
bcrank_p40D <- barcodeRanks(counts(p40D_sce))
bcrank_p40R <- barcodeRanks(counts(p40R_sce))
```

```{r plots}
plot_waterfall(bcrank_p16R,"p16R_sce")
plot_waterfall(bcrank_p16D,"p16D_sce")
plot_waterfall(bcrank_p40D,"p40D_sce")
plot_waterfall(bcrank_p40R,"p40R_sce")
plot_diag(p16D_sce)
plot_diag(p16R_sce)
plot_diag(p40D_sce)
plot_diag(p40R_sce)

p16D_sce <- empty_filt(p16D_sce) #dim: 5446 4260 
p16R_sce <- empty_filt(p16R_sce) #dim: 5408 5238
p40D_sce <- empty_filt(p40D_sce) #dim: 5363 3656
p40R_sce <- empty_filt(p40R_sce) #dim: 5405 3943 

```

### Percell QC

```{r}

p16D_sce <- add_qc(p16D_sce) 
p16R_sce <- add_qc(p16R_sce) 
p40D_sce <- add_qc(p40D_sce)  
p40R_sce <- add_qc(p40R_sce) 

plot_hist(p16D_sce)
plot_hist(p16R_sce)
plot_hist(p40D_sce)
plot_hist(p40R_sce)

saveRDS(p16D_sce,"p16D_empty_rm_QC.rds")
saveRDS(p16R_sce,"p16R_empty_rm_QC.rds")
saveRDS(p40D_sce,"p40D_empty_rm_QC.rds")
saveRDS(p40R_sce,"p40R_empty_rm_QC.rds")
```

### Normalising for cell Biases

This step adds the lognormalised counts matrix to sce object.

```{r}
p16D_sce2 <- normalise(p16D_sce,"p16D_sce2")
p16R_sce2 <- normalise(p16R_sce,"p16R_sce2")
p40D_sce2 <- normalise(p40D_sce,"p40D_sce2")
p40R_sce2 <- normalise(p40R_sce,"p40R_sce2")

plot(p16D_sce2_dec$mean, p16D_sce2_dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(p16D_sce2_dec)$trend(x), col="blue", add=TRUE)
plot(p16R_sce2_dec$mean, p16R_sce2_dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(p16R_sce2_dec)$trend(x), col="blue", add=TRUE)
plot(p40D_sce2_dec$mean, p40D_sce2_dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(p40D_sce2_dec)$trend(x), col="blue", add=TRUE)
plot(p40R_sce2_dec$mean, p40R_sce2_dec$total, xlab="Mean log-expression", ylab="Variance")
curve(metadata(p40R_sce2_dec)$trend(x), col="blue", add=TRUE)


p16D.top <- getTopHVGs(p16D_sce2_dec, prop=0.1)
p16R.top <- getTopHVGs(p16R_sce2_dec, prop=0.1)
p40D.top <- getTopHVGs(p40D_sce2_dec, prop=0.1)
p40R.top <- getTopHVGs(p40R_sce2_dec, prop=0.1)

p16D_sce2 <- dim_reduction(p16D_sce2,p16D.top, p16D_sce2_dec)
p16R_sce2 <- dim_reduction(p16R_sce2,p16R.top, p16R_sce2_dec)
p40D_sce2 <- dim_reduction(p40D_sce2,p40D.top, p40D_sce2_dec)
p40R_sce2 <- dim_reduction(p40R_sce2,p40R.top, p40R_sce2_dec)

# 
#plotUMAP(sce2, colour_by = "label")
#plotTSNE(sce2, colour_by = "label")
#top.dec <- subset(dec, rownames(dec) %in% top.hvgs)
#plotExpression(sce2, features=rownames(top.dec)[1:10])
```

### Remove doublets

```{r}

p16D_clean <-remove_doublet(p16D_sce2,p16D.top)
p16R_clean <-remove_doublet(p16R_sce2,p16R.top)
p40D_clean <-remove_doublet(p40D_sce2,p40D.top)
p40R_clean <-remove_doublet(p40R_sce2,p40R.top)


saveRDS(p16D_clean, "p16D_clean.rds")
saveRDS(p16R_clean, "p16R_clean.rds")
saveRDS(p40D_clean, "p40D_clean.rds")
saveRDS(p40R_clean, "p40R_clean.rds")
```

## Seurat

```{r}
rownames(p16D_clean) <- gsub("pf3d7_","",rownames(p16D_clean))
rownames(p16R_clean) <- gsub("pf3d7_","",rownames(p16R_clean))
rownames(p40D_clean) <- gsub("pf3d7_","",rownames(p40D_clean))
rownames(p40R_clean) <- gsub("pf3d7_","",rownames(p40R_clean))
p16D_clean$Sample <- "16D"
p16R_clean$Sample <- "16R"
p40D_clean$Sample <- "40D"
p40R_clean$Sample <- "40R"

p16D_clean$batch <- "16hpi"
p16R_clean$batch <- "16hpi"
p40D_clean$batch <- "40hpi"
p40R_clean$batch <- "40hpi"

p16D.seurat <- as.Seurat(p16D_clean)
p16R.seurat <- as.Seurat(p16R_clean)
p40D.seurat <- as.Seurat(p40D_clean)
p40R.seurat <- as.Seurat(p40R_clean)

#loading Malaria data
dat.m = read.csv("~/Desktop/Plasmodium_scRNAseq/malaria cell atlas data/pf10xIDC_counts.csv", header = TRUE ,row.names = 1)

mca.seurat <- CreateSeuratObject(counts = dat.m, project = "Malaria-Cell-Atlas")
mca.seurat[["percent.mt"]] <- PercentageFeatureSet(object = mca.seurat, pattern = "mal")
mca.seurat$Sample <- "MCA"
mca.seurat$batch <- "MCA"


saveRDS(sample.list,"sample.list.rds")
```
The `as.seurat` function do not automatically adds the `nCount_RNA` and the `nFeature_RNA` columns. Here is a quick trick for that: [Source](https://www.biostars.org/p/445065/).
Note: I did this after the entire analysis but this must be done here. Since the QC step is performed on raw counts we must change the default assay to RNA

```{r}
#sce <- readRDS("plasmodium.combined.with.mca.rds")
GetAssayData(sce, slot = "counts") ## This won't give you any matrix because default assay is integrated
DefaultAssay(seurat.obj) <- "RNA"
GetAssayData(seurat.obj,slot = "counts")
colSums(GetAssayData(seurat.obj,slot = "counts"))
seurat.obj$nCount_RNA <- colSums(GetAssayData(seurat.obj,slot = "counts"))
seurat.obj$nFeature_RNA <- colSums(GetAssayData(seurat.obj,slot = "counts")>0)
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "mal")
temp <- seurat.obj[,seurat.obj@meta.data$Sample=="16D"|seurat.obj@meta.data$Sample=="16R"|seurat.obj@meta.data$Sample=="40D"|seurat.obj@meta.data$Sample=="40R"]
VlnPlot(temp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "Sample")
p1 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Sample")
p2 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "Sample")
p1+p2

```

### Without SCTransform
```{r}
sample.list <- list(p16D=p16D.seurat,p16R=p16R.seurat,p40D=p40D.seurat,p40R=p40R.seurat,mca=mca.seurat)
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
#plasmodium.combined <- ScaleData(plasmodium.combined, verbose = FALSE)
plasmodium.combined <- ScaleData(plasmodium.combined, verbose = FALSE)
plasmodium.combined <- RunPCA(plasmodium.combined, npcs = 30, verbose = FALSE)

#p.combined.harmony <- RunHarmony(plasmodium.combined,assay.use="integrated", c("Sample", "batch"),plot_convergence = TRUE,verbose = TRUE)

p.combined.harmony <- RunHarmony(plasmodium.combined,assay.use="integrated","Sample",plot_convergence = TRUE,verbose = TRUE)
#harmony_embeddings <- Embeddings(p.combined.harmony, 'harmony')
p1 <- VlnPlot(object = plasmodium.combined, features = "PC_1", group.by = "Sample", pt.size = .001)
p2 <- VlnPlot(object = p.combined.harmony, features = "harmony_1", group.by = "Sample",pt.size = .001)
cowplot::plot_grid(p1,p2)

p.combined.harmony <- p.combined.harmony%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

#plasmodium.combined <- RunUMAP(plasmodium.combined, reduction = "pca", dims = 1:30)
#plasmodium.combined <- FindNeighbors(plasmodium.combined, reduction = "pca", dims = 1:30)
#plasmodium.combined <- FindClusters(plasmodium.combined, resolution = 0.5)

p1 <- DimPlot(p.combined.harmony, reduction = "umap", group.by = "Sample", cols = c('16D' = '#72147e', '16R' = '#fa9905', '40D' = '#04009a', '40R' = '#ff5200', 'MCA'='grey'))
p2 <- DimPlot(p.combined.harmony, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(p.combined.harmony, reduction = "umap", group.by = "batch", cols = c('16hpi' = '#72147e', '40hpi' = '#fa9905', 'MCA'='grey'))
p1+p2+p3
#saveRDS(p.combined.harmony,"plasmodium.combined.harmony.with.mca.rds")
saveRDS(p.combined.harmony,"plasmodium.combined.with.mac.harmony.sample_covariate.rds")
```

### With SCTransform

```{r}
sample.list_sct <- list(p16D=p16D.seurat,p16R=p16R.seurat,p40D=p40D.seurat,p40R=p40R.seurat,mca=mca.seurat)

#Run this step outsiide rmd file
sample.list_sct <- lapply(X = sample.list_sct, FUN =SCTransform)
features_sct <- SelectIntegrationFeatures(object.list = sample.list_sct,nfeatures = 2000)
sample.list_sct <- PrepSCTIntegration(object.list = sample.list_sct, anchor.features = features_sct)
plasmodium.anchors_sct <- FindIntegrationAnchors(object.list = sample.list_sct,normalization.method = "SCT", anchor.features = features_sct)
plasmodium.combined_sct <- IntegrateData(anchorset = plasmodium.anchors_sct,normalization.method = "SCT")
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
#DefaultAssay(plasmodium.combined_sct) <- "integrated"

# Run the standard workflow for visualization and clustering
plasmodium.combined_sct <- RunPCA(plasmodium.combined_sct, npcs = 30, verbose = FALSE)
#p.combined.sct.harmony <- RunHarmony(plasmodium.combined_sct, assay.use="integrated", c("Sample", "batch"),plot_convergence = TRUE,verbose = TRUE)
p.combined.sct.harmony <- RunHarmony(plasmodium.combined_sct, assay.use="integrated", "Sample",plot_convergence = TRUE,verbose = TRUE)

#harmony_embeddings_sct <- Embeddings(p.combined.sct.harmony, 'harmony')

p1 <- VlnPlot(object = p.combined.sct.harmony, features = "PC_1", group.by = "Sample", pt.size = .001)
p2 <- VlnPlot(object = p.combined.sct.harmony, features = "harmony_1", group.by = "Sample",pt.size = .001)
cowplot::plot_grid(p1,p2)

p.combined.sct.harmony <- p.combined.sct.harmony%>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = 0.5)

#plasmodium.combined_sct <- RunUMAP(plasmodium.combined_sct, reduction = "pca", dims = 1:30)
#plasmodium.combined_sct <- FindNeighbors(plasmodium.combined_sct, reduction = "pca", dims = 1:30)
#plasmodium.combined_sct <- FindClusters(plasmodium.combined_sct, resolution = 0.5)

p1 <- DimPlot(p.combined.sct.harmony, reduction = "umap", group.by = "Sample",cols = c('16D' = '#72147e', '16R' = '#fa9905', '40D' = '#04009a', '40R' = '#ff5200', 'MCA'='grey'))
p2 <- DimPlot(p.combined.sct.harmony, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(p.combined.sct.harmony, reduction = "umap", group.by = "batch",cols = c('16hpi' = '#72147e', '40hpi' = '#fa9905', 'MCA'='grey'))
p1+p2+p3
#saveRDS(p.combined.sct.harmony,"plasmodium.combined.with.mca.harmony.sct.rds")
saveRDS(p.combined.sct.harmony,"plasmodium.combined.with.mac.sct.harmony.sample_covariate.rds")
```
### Cell type annotation of seurat object

```{r}

mca.pheno = read.csv("/home/pichkari/Desktop/Plasmodium_scRNAseq/malaria cell atlas data/pf10xIDC_pheno.csv")

#Since we are adding metadata about cell types of MCA after running seurat, the barcodes are appended with sample group number (here _5 for MCA)
mca.pheno$ids <- paste0(mca.pheno$X,"_5")


p <- readRDS("plasmodium.combined.with.mac.harmony.batch_covariate.rds")
p <- map(p)
saveRDS(p,"plasmodium.combined.with.mac.harmony.batch_covariate.rds")
p <- readRDS("plasmodium.combined.with.mac.sct.harmony.batch_covariate.rds")
p <- map(p)
saveRDS(p,"plasmodium.combined.with.mac.sct.harmony.batch_covariate.rds")

p <- readRDS("./RDS_files/plasmodium.combined.with.mca.rds")
p <- map(p)
saveRDS(p,"./RDS_files/plasmodium.combined.with.mca.rds")

```


### Identify conserved cell type markers

```{r}
library(metap)

# We will use the sct transformed harmony seurat object and non-sct-non-harmony seurat object
p.combined.with.mca <- readRDS("./RDS_files/plasmodium.combined.with.mca.rds")
p1 <- DimPlot(p.combined.with.mca, reduction = "umap", group.by = "Sample",cols = c('16D' = '#72147e', '16R' = 'green', '40D' = '#04009a', '40R' = '#ff5200', 'MCA'='grey'))
p2 <- DimPlot(p.combined.with.mca, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(p.combined.with.mca, reduction = "umap", group.by = "batch",cols = c('16hpi' = '#72147e', '40hpi' = '#fa9905', 'MCA'='grey'))
p4 <- DimPlot(p.combined.with.mca, reduction = "umap", label = TRUE, repel = TRUE,group.by = "bulk")
library(ggpubr)
ggarrange(p1,p3,p2,p4,labels = c("A","B","C","D"),nrow = 2, ncol = 2)

#############################################
p.combined.with.mca.sct.harmomy.batch <- readRDS("plasmodium.combined.with.mac.sct.harmony.batch_covariate.rds")
p1 <- DimPlot(p.combined.with.mca.sct.harmomy.batch, reduction = "umap", group.by = "Sample",cols = c('16D' = '#72147e', '16R' = 'green', '40D' = '#04009a', '40R' = '#ff5200', 'MCA'='grey'))
p2 <- DimPlot(p.combined.with.mca.sct.harmomy.batch, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(p.combined.with.mca.sct.harmomy.batch, reduction = "umap", group.by = "batch",cols = c('16hpi' = '#72147e', '40hpi' = '#fa9905', 'MCA'='grey'))
p4 <- DimPlot(p.combined.with.mca.sct.harmomy.batch, reduction = "umap", label = TRUE, repel = TRUE,group.by = "bulk")
library(ggpubr)
ggarrange(p1,p3,p2,p4,labels = c("A","B","C","D"),nrow = 2, ncol = 2)
```

### Plotting percentage of cells per cluster per batch
```{r}
## Usage
## On seurat obj Without sct and harmony
t <- proportion_calc(sobj = p.combined.with.mca,rowsum = TRUE,rows = "seurat_clusters", cols = "Sample")
proportion_plot(t,xlab = "Clusters",ylab = "Percentage of cells.per.sample.per.cluster", legend.title = "Sample")
## On seurat object with sct and harmony(corrected for batch only)
t <- proportion_calc(sobj = p.combined.with.mca.sct.harmomy.batch,rowsum = TRUE,rows = "seurat_clusters", cols = "Sample")
proportion_plot(t,xlab = "Clusters",ylab = "Percentage of cells.per.sample.per.cluster", legend.title = "Sample")

##Plot mca cell type percentageper cluster
t <- proportion_calc(sobj = p.combined.with.mca,rowsum = TRUE,rows = "seurat_clusters", cols = "bulk",retain.mca = TRUE)
proportion_plot(t,xlab = "Clusters",ylab = "Percentage of cells.per.celltype.per.cluster", legend.title = "Cell Type")
```
On the basis of percent contribution of cell types within each cluster, cluster were annotated. Cluster belonged to that cell type which contribute more percentage of mca cells. 
```{r}
# annotating clusters
t <- RenameIdents(p.combined.with.mca, `0` = "Ring", `1` = "Early_Troph", `2` = "Early_Troph",`3` = "Ring", `4` = "Late_Troph", `5` = "Late_Troph", `6` = "Late_Troph", `7` = "Late_Troph", `8` = "Schizont", `9` = "Schizont",`10` = "Late_Troph", `11` = "Late_Troph", `12` = "Schizont", `13` = "Ring", `14` = "Late_Troph")

## The new cluster annotations are stored in active.ident slot of seurat. So to add it to metadata table as new column, we do: 
DimPlot(t, label = TRUE)
#identical(rownames(t@meta.data),rownames(data.frame(t@active.ident)))
## Saving these cluster classification in annot column
p.combined.with.mca@meta.data$annot <- data.frame(t@active.ident)$t.active.ident

t <- proportion_calc(sobj = p.combined.with.mca,rowsum = FALSE,rows = "annot", cols = "Sample",retain.mca = FALSE)
proportion_plot(t(t),xlab = "Samples",ylab = "Percentage of cells.per.celltype", legend.title = "Cell Type")
######### Classifying the SCT harmonized object
t=RenameIdents(p.combined.with.mca.sct.harmomy.batch, `0` = "Ring", `1` = "Late_Troph", `2` = "Early_Troph", 
    `3` = "Late_Troph", `4` = "Late_Troph", `5` = "Early_Troph", `6` = "Schizont", `7` = "Schizont", `8` = "Late_Troph", `9` = "Late_Troph", 
    `10` = "Late_Troph", `11` = "Schizont", `12` = "Ring", `13` = "Late_Troph")
DimPlot(t, label = TRUE)
p.combined.with.mca.sct.harmomy.batch@meta.data$annot <- data.frame(t@active.ident)$t.active.ident

t <- proportion_calc(sobj = p.combined.with.mca.sct.harmomy.batch,rowsum = FALSE,rows = "annot", cols = "Sample",retain.mca = FALSE)
proportion_plot(t(t),xlab = "Samples",ylab = "Percentage of cells.per.celltype", legend.title = "Cell Type")
```
### Var gene expression per cell
```{r}
#ensembl_var_genes <- read.csv("~/rohit/kaust_proj/var_gene.csv", header = F)
kaust_var_genes <- readr::read_tsv("~/rohit/kaust_proj/kaust_vargene_list.tsv", col_names = TRUE)

## Dont use subset function of seurat


## The function  below spits out seurat object with extra column in metadata named gene_exp (var_gene expression for provided gene list)
#eg:percentage expression per cell of genes in geneList
temp <- gene_exp_per(p.combined.with.mca,kaust_var_genes$`Gene ID`,"gene_expr")

### Function for Number of genes in geneList expressed per cell


```
### Finding Conserved Markers

```{r}
p.combined.with.mca <- readRDS("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/p.combined.with.mca.app.rds")
p.combined.with.mca$sample.cluster.idents <- paste(Idents(p.combined.with.mca),p.combined.with.mca$Sample,sep = "_")

## Remove MCA
p.combined.with.mca <-p.combined.with.mca[,p.combined.with.mca@meta.data$Sample=="16D"|p.combined.with.mca@meta.data$Sample=="16R"|p.combined.with.mca@meta.data$Sample=="40D"|p.combined.with.mca@meta.data$Sample=="40R"]

DefaultAssay(p.combined.with.mca) <- "RNA"

## Find Conserved Markers. To  calculate the genes that are conserved markers irrespective of condition (provided by column sample) in all cluster 0:14

conserved.markers <- function(sobj, clusters,group.by,min.pct = 0.50){
  for (i in clusters){
    markers <- Seurat::FindConservedMarkers(sobj, ident.1 = i, grouping.var = group.by, verbose = TRUE,min.pct = min.pct)
    write.csv(markers,file = paste0("cluster.",i,".conserved.markers.csv"))
    print(paste0("Writing the file ","cluster.",i,".conserved.markers.csv"))
    }
}
clusters <- c(0,1,2,3,4,6,7,10,11,13,14)
clusters_16hpi <- c(0,1,2,3,4,6,7,10,11,13,14)
clusters_40hpi <- c(0,1,2,3,4,5,6,7,8,9,11,12,13,14)
conserved.markers(sobj = p.combined.with.mca,clusters = clusters,group.by = "Sample",min.pct = 0.50)
## Finding DE Genes across conditions
de.sobj <- p.combined.with.mca

Idents(de.sobj) <- de.sobj$sample.cluster.idents

de.genes <- function(sobj,ident1,ident2,clusters,min.pct = 0.50, batch.name){
  for (i in clusters){
    mk <- FindMarkers(sobj,ident.1 = paste(i,ident1,sep = "_"), ident.2 =paste(i,ident2,sep = "_"), min.pct = min.pct)
    write.csv(mk,file = paste0("de_",i,"_",batch.name,".csv"))
    print(paste0("Writing the file ","de_",i,"_",batch.name,".csv"))
    }
}

de.genes(sobj = de.sobj,ident1 = "16D",ident2 = "16R",clusters = clusters_16hpi, min.pct = 0.50,batch.name = "16hpi")
de.genes(sobj = de.sobj,ident1 = "40D",ident2 = "40R",clusters = clusters_40hpi, min.pct = 0.50,batch.name = "40hpi")
```
# GO Term analysis
[source1](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/)
[source2](https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/)
```{r}
## Keep checking the clusterprofiler document because arguments keeps on changing
library(clusterProfiler)
library(org.Pf.plasmo.db)
library(enrichplot)
library(wordcloud)
pf.go <- org.Pf.plasmo.db
cluster0 <- read.csv("de_10_16D.csv")
# gene_list <- cluster0$avg_log2FC
# names(gene_list) <- gsub("-","_",cluster0$X)
# gene_list<-na.omit(gene_list)
# gene_list = sort(gene_list, decreasing = TRUE)
genes <- gsub("-","_",cluster0$X)
go <- enrichGO(gene=genes, 
             ont ="ALL",keyType="SYMBOL",
             OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
             qvalueCutoff = 0.10, 
             pAdjustMethod = "BH")
dotplot(go, showCategory=30)
barplot(go, 
        drop = TRUE, 
        showCategory = 30, 
        title = "GO Biological Pathways",
        font.size = 8)
#View(as.data.frame(go))
## Frequency
wcdf<-read.table(text=go$GeneRatio, sep = "/")[1]
## function terms
wcdf$term<-go[,3]
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
```
## Chec percentage of cells with exon2 silenced


```{r}
setwd("~/rohit/kaust_proj")
p16D <- readr::read_tsv("/home/pichkari/rohit/kaust_proj/16D_with_exon2_reads.sam",col_names = FALSE)
p16R <- readr::read_tsv("/home/pichkari/rohit/kaust_proj/16R_with_exon2_reads.sam",col_names = FALSE)
p40D <- readr::read_tsv("/home/pichkari/rohit/kaust_proj/40D_with_exon2_reads.sam",col_names = FALSE)
p40R <- readr::read_tsv("/home/pichkari/rohit/kaust_proj/40R_with_exon2_reads.sam",col_names = FALSE)

filter <- function(df,nh=12,cb,append){
  all <- df[,c(nh,cb)]
  all <- unique(all)
  all <- all[all[1]=="NH:i:1",]
  cb <- list()
for (i in 2:ncol(all)){
  d <- all[,i]
  print(i)
  colnames(d) <- "col"
  cb[i] <- all[grep("CB:Z:",d$col),][i]
}
  cb_all <- data.frame("cb"=unlist(cb))[1]
  cb_all <- gsub("CB:Z:","",cb_all$cb)
  cb_all <- gsub("-1",paste0("-1",append),cb_all)
  return(cb_all)
  }

p16D.ids <- filter(p16D,nh=12,cb=c(21,22,23),append = "_1")

p16R.ids <- filter(p16R,nh=12,cb=c(21,22,23),append = "_2")

p40D.ids <- filter(p40D,nh=12,cb=c(21,22,23,25),append = "_3")

p40R.ids <- filter(p40R,nh=12,cb=c(21,22,23,25),append = "_4")

p.combined.with.mca <- readRDS("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/p.combined.with.mca.app.rds")
p.combined.with.mca$sample.cluster.idents <- paste(Idents(p.combined.with.mca),p.combined.with.mca$Sample,sep = "_")

## Remove MCA
p.combined.with.mca <-p.combined.with.mca[,p.combined.with.mca@meta.data$Sample=="16D"|p.combined.with.mca@meta.data$Sample=="16R"|p.combined.with.mca@meta.data$Sample=="40D"|p.combined.with.mca@meta.data$Sample=="40R"]

DefaultAssay(p.combined.with.mca) <- "RNA"
p.subset <- p.combined.with.mca["PF3D7-1107800",]
p.subset$pfap2_exp_status <- "pfAP2-MRP Not Expressed"
pfap2_exp_cb <-names(p.subset@assays$RNA@counts[,colSums(p.subset@assays$RNA@counts)>0])

p.subset@meta.data[which(colnames(p.subset) %in% pfap2_exp_cb),]$pfap2_exp_status <- "pfAP2-MRP Expressed"
t <- table(p.subset$pfap2_exp_status,p.subset$sample.cluster.idents)

percentage <- function(t){
  t <- data.frame(t)
  v <- as.vector(unique(t$Var2))
  t$total <- "total"
  for(i in v){
  t[t$Var2==i,]$total <- colSums(t[t$Var2==i,][3])
  }
  t$percent <- (as.numeric(t$Freq))*100/as.numeric(t$total)
return(t)
}

t <- percentage(t)
subset.16R <- t[grep("16R",t$Var2),]
subset.16D <- t[grep("16D",t$Var2),]
ggplot(subset.16R, aes(fill=Var1, y=percent, x=Var2))+geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_viridis(discrete = T)+xlab("Clusters")+ylab("Percentage of 16R Cells in cluster expressing PfAP2-MRP")+size+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(fill="pfAp2-MRP")

p.subset$knockout_status <- "Exon2 Not Expressed"

p.subset@meta.data[which(colnames(p.subset) %in% p16D.ids),]$knockout_status <- "Exon2 Expressed"
p.subset@meta.data[which(colnames(p.subset) %in% p16R.ids),]$knockout_status <- "Exon2 Expressed"
p.subset@meta.data[which(colnames(p.subset) %in% p40D.ids),]$knockout_status <- "Exon2 Expressed"
p.subset@meta.data[which(colnames(p.subset) %in% p40R.ids),]$knockout_status <- "Exon2 Expressed"
#p.cells.expressing@meta.data[which(colnames(p.cells.expressing) %in% p40D.ids),]$knock_out_status <- "Exon2 Reads Present"

#p.cells.expressing@meta.data[which(colnames(p.cells.expressing) %in% p16R.ids),]$knock_out_status <- "Exon2 Reads Present"

#p.cells.expressing@meta.data[which(colnames(p.cells.expressing) %in% p16D.ids),]$knock_out_status <- "Exon2 Reads Present"#

t <- table(p.combined.with.mca$knock_out_status,p.combined.with.mca$seurat_clusters)

## Now that we have all cells that express AP2-MRP gene, we wish to see which cluster has maximum percentage of cells expressing this gene
t <- table(p.cells.expressing$Sample,p.cells.expressing$seurat_clusters)
t <- percentage(t)
ggplot(t, aes(fill=Var1, y=percent, x=Var2)) + 
    geom_bar(position="stack", stat="identity")+theme_bw()+scale_fill_viridis(discrete = T)+xlab("Sample")+ylab("No. of Cells Failed Knockout")+labs(fill="Status")+size

t <- table(p.combined.with.mca$knock_out_status,p.combined.with.mca$Sample)
t <- percentage(t)
t <- table(p.combined.with.mca$knock_out_status,p.combined.with.mca$sample.cluster.idents)
t <- percentage(t)
```


# Harmony with batch only Done
# MCA bulk mapping -> table() -> barplot by row Done
# Find clusters -> table -> percentage of each sample per cluster DOne
# Find conserved markers with batch only










