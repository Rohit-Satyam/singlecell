## Finding the Var genes:

We downloaded the gff file from PlasmoDB v37

```bash
grep PfEMP1 PlasmoDB-37_Pfalciparum3D7.gff | awk '{print $9}' | sort -u | cut -f '2' -d'=' | cut -f '1' -d';' | cut -f '1' -d'.' | sort -u > var_genes.csv
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
sobj.m$Sample <- "mca"
sobj.m$batch <- "MCA"

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
saveRDS(plasmodium.combined,"plasmodium.combined.with.mca.rds")
var_genes <- read.csv("~/rohit/kaust_proj/var_gene.csv")

# Finding common features
r2 <- var_genes$V1
r <- rownames(plasmodium.combined)
common <- intersect(r,r2)

#intersect(r,r2)
 #[1] "PF3D7-1000100" "PF3D7-0632500" "PF3D7-0412400" "PF3D7-1240600" "PF3D7-0412700" "PF3D7-1240400"
 #[7] "PF3D7-0712900" "PF3D7-1200600" "PF3D7-0712300" "PF3D7-0712000" "PF3D7-0421300" "PF3D7-0632800"
#[13] "PF3D7-1240900" "PF3D7-0712800" "PF3D7-0800300" "PF3D7-0426000" "PF3D7-0900100" "PF3D7-0712600"
#[19] "PF3D7-0223500" "PF3D7-0712400" "PF3D7-0420700" "PF3D7-1150400" "PF3D7-0711700" "PF3D7-1100200"
#[25] "PF3D7-0600200" "PF3D7-1041300" "PF3D7-1219300" "PF3D7-0800200" "PF3D7-0937800" "PF3D7-1200400"
#[31] "PF3D7-1240300" "PF3D7-0425800" "PF3D7-0421100" "PF3D7-0808700" "PF3D7-0808600" "PF3D7-1100100"
#[37] "PF3D7-1373500" "PF3D7-0733000" "PF3D7-0300100" "PF3D7-0617400" "PF3D7-0115700" "PF3D7-1200100"
#[43] "PF3D7-0809100" "PF3D7-0324900"

FeaturePlot(object = plasmodium.combined, features = common[1:12],min.cutoff = "q9" )
FeaturePlot(object = plasmodium.combined, features = common[13:24],min.cutoff = "q9" )
FeaturePlot(object = plasmodium.combined, features = common[25:40],min.cutoff = "q9" )
FeaturePlot(object = plasmodium.combined, features = common[41:47],min.cutoff = "q9" )

t <- plasmodium.combined@assays$RNA[common,] ## 44 Var Genes

t2 <- colSums(t)
t2 <- Matrix::colSums(t)
table(t2>0)

#FALSE  TRUE 
#11013 13613 

#11013 cells that do not express any of 44 var genes.

g <- t[,t2>0]
m <- colSums(g)
table(m>1)

#FALSE  TRUE 
 #88 13525 

 var_mat = (plasmodium.combined@assays$RNA)[common,]>0
 
 dim(var_mat)
 
 #dim(var_mat)
#[1]    44 24626

table(colSums(var_mat))

#    0     1     2     3     4     5     6     7     8     9    10    11    12    14 
#11013 10189  2326   626   250   115    58    21    13     5     4     4     1     1 

#vector with cells containing more than 1 var gene
vv = c( sum( (colSums(var_mat)>1)[which(plasmodium.combined$Sample=="16R")]),
sum( (colSums(var_mat)>1)[which(plasmodium.combined$Sample=="16D")]),
sum( (colSums(var_mat)>1)[which(plasmodium.combined$Sample=="40R")]),
sum( (colSums(var_mat)>1)[which(plasmodium.combined$Sample=="40D")]) )

vv


boxplot(
  colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="16R")],
  colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="16D")],
  colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="40R")],
  colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="40D")]
)

# proportion of cells expressing more than 1 var gene increases with hpi and in mutation

vv_var0=
  cbind(
    table((colSums(var_mat)>0)[which(plasmodium.combined$Sample=="16R")]),
    table((colSums(var_mat)>0)[which(plasmodium.combined$Sample=="16D")]),
    table((colSums(var_mat)>0)[which(plasmodium.combined$Sample=="40R")]),
    table((colSums(var_mat)>0)[which(plasmodium.combined$Sample=="40D")])
  )

for(i in 1:ncol(vv_var0))
  vv_var0[,i]=vv_var0[,i]/sum(vv_var0[,i])

barplot(vv_var0,names.arg = c("16R","16D","40R","40D"),ylab="% cells expressing var genes")

vv_var1more=
  cbind(
  table(colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="16R")]>1),
  table(colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="16D")]>1),
  table(colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="40R")]>1),
  table(colSums(var_mat)[which( colSums(var_mat)>0 & plasmodium.combined$Sample=="40D")]>1)
)

for(i in 1:ncol(vv_var1more))
  vv_var1more[,i]=vv_var1more[,i]/sum(vv_var1more[,i])

barplot(vv_var1more,names.arg = c("16R","16D","40R","40D"),ylab="% cells expressing 1+ var genes")

#fisher test on cells with more than 1 var gene in the context of R,D and 16,40
fisher.test(matrix(vv,nrow = 2,byrow = T))
# suggests that with later hpi and mutation the #var genes and #cells expressing var genes significantly increases


# Conclusion with later hpi lesser proportion of cells express var genes
# However, at later hpi more number of var genes per cells are expressed
# At 40 hpi D, both the number of cells expressing var genes and cells expressing more than 1 var genes is greater

## Make Heatmap of the var genes from KAUST
readr::read_tsv("~/rohit/kaust_proj/kaust_vargene_list.tsv", col_names = TRUE)
r3 <- kaust_var_genes$`Gene ID`
common_kaust_var_genes <- intersect(r,r3)
DoHeatmap(subset_sobj, features = common_kaust_var_genes[1:5],  group.by = "Sample")+theme(axis.text.y = element_text(color="black", size=20))

```


