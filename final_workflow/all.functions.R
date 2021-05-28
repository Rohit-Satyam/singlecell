suppressPackageStartupMessages({
  # library(scran) # giving error
  library(viridis)
  library(multtest)
  #library(harmony)
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

### Filtering for non-expressed genes
filter <- function(sobj){
  keep_features <- rowSums(counts(sobj) > 0) > 0
  sobj <- sobj[keep_features, ]
  name <- rownames(sobj)
  hg <- grep("hg38", rownames(sobj))
  name <- name[-hg]
  sobj <-sobj[name,]
  return(sobj)
}
### Diagnosing for empty cells (Columns)
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

### Percell QC
add_qc <- function(sce){
  is.mito <- grep("mal", rownames(sce))
  sce <-  addPerCellQC(sce, subsets=list(Mito=is.mito))
  return(sce)}

plot_hist <- function(sce){
  hist(sce$subsets_Mito_percent, breaks=20, col="grey80",
       xlab="Proportion of reads in mitochondrial genes")
}
#Here dt is seurat object derieved df with Subset_Mito_percent and Sample column
## p3 <- ggplot(dt[which(dt$sample=="16D"),], aes(x=mito_percent, color=sample, fill = sample )) + 
#geom_histogram(alpha=0.5, position="identity")+theme_bw()+theme(axis.text.x = element_text(face="bold", color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=12),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.title = element_text(face="bold",color = "black", size = 12),legend.text = element_text(face = "bold",color = "black",size = 12))+xlab("Per Cell mitochondrial read percentage")

### Normalising for cell Biases

normalise <- function(sobj, name){
  clusters <- quickCluster(sobj, BSPARAM=IrlbaParam())
  print(table(clusters))
  sce <- computeSumFactors(sobj, min.mean=0.1, cluster=clusters)
  summary(sizeFactors(sce))
  sce <- logNormCounts(sce)
  assign(paste0(name,"_dec"), modelGeneVarByPoisson(sce), envir = .GlobalEnv)
  return(sce)
}

dim_reduction <- function(sobj, top, dec){
  # Evaluate PCs
  sce2 <- denoisePCA(sobj, subset.row = top, technical = dec, BSPARAM=IrlbaParam())
  # make TSNE plot
  sce2 <- runTSNE(sce2, dimred = "PCA")
  # make UMAP plot
  sce2 <- runUMAP(sce2, dimred = "PCA")
  g <- buildSNNGraph(sce2, k = 10, use.dimred = "PCA")
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(sce2) <- factor(clust)
  return(sce2)
}

### Remove Doublets
remove_doublet <- function(sce2,top){
  dbl.dens <- scDblFinder::computeDoubletDensity(sce2,  d=ncol(reducedDim(sce2)),subset.row=top)
  sce2$DoubletScore <- dbl.dens
  cut_off <- quantile(sce2$DoubletScore, 0.95)
  sce2$isDoublet <- c("no", "yes")[factor(as.integer(sce2$DoubletScore >= cut_off),levels = c(0, 1))]
  sce_clean <- sce2[,sce2$isDoublet == "no"]
  return(sce_clean)
}

### Map MCA data on Seurat object
map <- function(sobj)
{
  sobj$rownams <- rownames(sobj@meta.data)
  t <- S4Vectors::merge(x=sobj@meta.data, y=mca.pheno[,c("ids","bulk")], by.x='rownams', by.y='ids', all.x=TRUE, sort=FALSE)%>% arrange(factor(rownams, levels = sobj@meta.data$rownams))%>%`rownames<-`(.$"rownams")
  print(identical(t$rownams,sobj@meta.data$rownams))
  sobj@meta.data <- t
  return(sobj)
}

## Function to compute proportion. Plotting percentage of cells per cluster per batch
proportion_calc <- function(sobj, rowsum=TRUE, rows, cols, retain.mca=FALSE ){
  t <-table(sobj@meta.data[,rows],sobj@meta.data[,cols])
  ifelse(retain.mca==TRUE,t,t<-t[,1:ncol(t)-1])
  if(rowsum==TRUE){
    for(i in 1:nrow(t)) {t[i,]=t[i,]*100/sum(t[i,])}
    
  }else{
    for(i in 1:ncol(t)) t[,i]=t[,i]*100/sum(t[,i])
  }
  
  return(t)
  
}

proportion_plot <- function(table,xlab,ylab,legend.title){
  ggplot(data.frame(table), aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat="identity")+scale_fill_viridis(discrete = T)+theme_bw()+theme(axis.text.x = element_text(face="bold", color="black", size=12),axis.text.y = element_text(face="bold", color="black", size=12),axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),legend.title = element_text(color = "black", size = 12),legend.text = element_text(color = "black",size = 12))+xlab(xlab)+ylab(ylab)+ labs(fill=legend.title)
}

#Subset the seurat object and remove MCAcells
# Code below was originally used:
# gene_exp_per <- function(sobj,genelist, colname){
#   q <- sobj[,sobj@meta.data$Sample=="16D"|sobj@meta.data$Sample=="16R"|sobj@meta.data$Sample=="40D"|sobj@meta.data$Sample=="40R"]
#   t <-data.frame(colSums(q@assays$RNA[which( rownames(q@assays$RNA) %in% genelist),])/colSums(q@assays$RNA))
#   colnames(t) <- colname
#   t$names <- rownames(t)
#   sobj@meta.data <-S4Vectors::merge(x=sobj@meta.data, y=t, by.x='rownams', by.y='names', all.x=TRUE, sort=FALSE)%>% arrange(factor(rownams, levels = sobj@meta.data$rownams))%>%`rownames<-`(.$"rownams")
#   return(sobj)
# }

## But for app we will use this one:
gene_exp_per <- function(sobj,genelist, colname){
  q <- sobj
  t <-data.frame(colSums(q@assays$RNA[which( rownames(q@assays$RNA) %in% genelist),])/colSums(q@assays$RNA))
  colnames(t) <- colname
  t$names <- rownames(t)
  sobj@meta.data <-S4Vectors::merge(x=sobj@meta.data, y=t, by.x='rownams', by.y='names', all.x=TRUE, sort=FALSE)%>% arrange(factor(rownams, levels = sobj@meta.data$rownams))%>%`rownames<-`(.$"rownams")
  return(sobj)
}




