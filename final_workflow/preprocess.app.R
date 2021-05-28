## Making app ready rds seurat object
seurat.obj <- readRDS("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/RDS_files/plasmodium.combined.with.mca.rds")
mca.pheno = read.csv("/home/pichkari/Desktop/Plasmodium_scRNAseq/malaria cell atlas data/pf10xIDC_pheno.csv")
mca.pheno$ids <- paste0(mca.pheno$X,"_5")
source('~/rohit/kaust_proj/sample_swap_corrected/all_functions.R')
kaust_var_genes <- readr::read_tsv("~/rohit/kaust_proj/kaust_vargene_list.tsv", col_names = TRUE)
surfins <- kaust_var_genes[grep("SURFIN",kaust_var_genes$`Protein Product`, ignore.case = TRUE),]
stevor <- kaust_var_genes[grep("stevor",kaust_var_genes$`Protein Product`, ignore.case = TRUE),]
rifin <- kaust_var_genes[grep("rifin",kaust_var_genes$`Protein Product`, ignore.case = TRUE),]
pfmc2TM<- kaust_var_genes[grep("Pfmc-2TM",kaust_var_genes$`Protein Product`,ignore.case = T),]
var_genes <- kaust_var_genes[grep("PfEMP1",kaust_var_genes$`Protein Product`,ignore.case = T),]

DefaultAssay(seurat.obj) <- "RNA"
seurat.obj$nCount_RNA <- colSums(GetAssayData(seurat.obj,slot = "counts"))
seurat.obj$nFeature_RNA <- colSums(GetAssayData(seurat.obj,slot = "counts")>0)
seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "mal")
DefaultAssay(seurat.obj) <- "integrated"
t <- RenameIdents(seurat.obj, `0` = "Ring", `1` = "Early_Troph", `2` = "Early_Troph",`3` = "Ring", `4` = "Late_Troph", `5` = "Late_Troph", `6` = "Late_Troph", `7` = "Late_Troph", `8` = "Schizont", `9` = "Schizont",`10` = "Late_Troph", `11` = "Late_Troph", `12` = "Schizont", `13` = "Ring", `14` = "Late_Troph")
seurat.obj@meta.data$annot <- data.frame(t@active.ident)$t.active.ident
########################### Processing for Var Gene Analysis ###
#var gene expression fraction per cell
seurat.obj <- gene_exp_per(seurat.obj,var_genes$`Gene ID`,"var_genes_percent")
seurat.obj <- gene_exp_per(seurat.obj,rifin$`Gene ID`,"rifin_percent")
seurat.obj <- gene_exp_per(seurat.obj,stevor$`Gene ID`,"stevor_percent")
seurat.obj <- gene_exp_per(seurat.obj,surfins$`Gene ID`,"surfin_percent")
seurat.obj <- gene_exp_per(seurat.obj,pfmc2TM$`Gene ID`,"pfc2tm_percent")

raw_mat = seurat.obj@assays$RNA@counts
raw_mat <- raw_mat[var_genes$`Gene ID`,]
seurat.obj@meta.data$var_num =colSums(raw_mat !=0)
seurat.obj@meta.data$var_num_per_tot_feature = seurat.obj@meta.data$var_num/seurat.obj@meta.data$nFeature_RNA
saveRDS(seurat.obj,"p.combined.with.mca.app.rds")
