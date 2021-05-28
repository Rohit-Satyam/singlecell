## Loading 16D count matrix
## Matrix of Var_genes from Cell Ranger with our python script run
mx.counts.py.run<- read.csv("/home/pichkari/rohit/kaust_proj/mapping_check/subset_10x-bams-Var-genes/16R_count_var.genes.csv", header = T)
mx.counts.py.run.new<- read.csv("/home/pichkari/rohit/kaust_proj/mapping_check/subset_10x-bams-Var-genes/16R_count_new.csv", header = T)

## Matrix of var genes filtered according to STAR Read IDs only and 
#80% reads mapping uniquely criteria
mx.counts.py.run.star.80 <- read.csv("/home/pichkari/rohit/kaust_proj/mapping_check/subset_10X-derieved-80-percent-Uniqely-mapped-reads/16R_count_80.csv", header = T)
## Matrix of var genes filtered according to STAR Read IDs only
mx.counts.py.run.star <- read.csv("/home/pichkari/rohit/kaust_proj/mapping_check/subset_10X-STAR-derieved-bams/16R_count.csv", header = T)
 
## Make matrix rownames and colnames comparable with that of Seurat
alter_matx <- function(mat_new){
mat_new$pseudoname <- gsub("_","-",mat_new$pseudoname)
mat_new$pseudoname <- gsub("pf3d7-","",mat_new$pseudoname)
rownames(mat_new) <- mat_new$pseudoname
mat_new <- mat_new[,-1]
colnames(mat_new) <- paste0(colnames(mat_new),"_2")
return(mat_new)
}

mx.counts.py.run <- alter_matx(mx.counts.py.run)
mx.counts.py.run.new <- alter_matx(mx.counts.py.run.new)
mx.counts.py.run.star <- alter_matx(mx.counts.py.run.star)
mx.counts.py.run.star.80 <- alter_matx(mx.counts.py.run.star.80)
sobj <- p.combined.with.mca[,p.combined.with.mca@meta.data$Sample=="16R"]
t <- data.frame(sobj@assays$RNA@counts)

## Compute the correlation between the seurat object t for 16D and given matrix but
## sampling 500 cells randomly from the colnames
correlation <- function(mat_new){
  mat_new <- mat_new[,intersect(colnames(mat_new),colnames(t))]
  mat_new <- mat_new[,sample(ncol(mat_new), size = 500), drop = FALSE]
common_row <- intersect(rownames(mat_new), rownames(t))
common_col <- intersect(colnames(mat_new), colnames(t))
mat.subset <- as.matrix(mat_new[common_row,common_col])
t.subset <- as.matrix(t[common_row,common_col])
## Column wise correlation
col.cor <- sapply(seq.int(dim(mat.subset)[2]), function(i) cor(mat.subset[,i],t.subset[,i]))
return(col.cor)
}

temp1 <- correlation(mx.counts.py.run)
hist(temp1, ylab="Number of Randomly Sampled Cells ",xlab="Corr. Coeff.", breaks = 100)
temp2 <- correlation(mx.counts.py.run.star)
hist(temp2, ylab="Number of Randomly Sampled Cells ",xlab="Corr. Coeff.", breaks = 100)
temp3 <- correlation(mx.counts.py.run.star.80)
hist(temp3, ylab="Number of Randomly Sampled Cells ",xlab="Corr. Coeff.", breaks = 100)


row.cor <- sapply(seq.int(dim(mat.subset)[1]), function(i) cor(mat.subset[i,],t.subset[i,]))

## Select cells expressing more thaan one Var genes


one.or.more <- function(mtx){
  return(colSums(mtx !=0))
}

w1 <- one.or.more(mx.counts.py.run)
table(w1)
t_var <- t[var_genes$`Gene ID`,]
w2 <- one.or.more(t_var)
table(w2)
# 
genes <- function(matx,cellbcode){
  return(matx %>% select(.,cellbcode) %>%filter(.[,1]> 0))
}
  
}

