library(shiny)
library(shinydashboard)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(metap)
library(viridis)
library(clusterProfiler)
library(org.Pf.plasmo.db)
library(enrichplot)
library(wordcloud)
set.seed(1000)
pf.go <- org.Pf.plasmo.db
clusters_16hpi <- c(0,1,2,3,4,6,7,11,13,14)
clusters_40hpi <- c(0,1,2,3,4,5,6,7,8,9,11,12,13,14)

size <- theme_bw()+theme(text = element_text(face="bold",color="black",size = 20),legend.title = element_text(face="bold",color = "black", size = 15),legend.text = element_text(face = "bold",color = "black",size = 15))
############ Change Path Here #############################

seurat.obj <- readRDS("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/p.combined.with.mca.app.rds")
seurat.obj$bulk <- tidyr::replace_na(seurat.obj$bulk,"Kaust Single Cells")
kaust_var_genes <- readr::read_tsv("~/rohit/kaust_proj/kaust_vargene_list.tsv", col_names = TRUE)
surfins <- kaust_var_genes[grep("SURFIN",kaust_var_genes$`Protein Product`, ignore.case = TRUE),]
stevor <- kaust_var_genes[grep("stevor",kaust_var_genes$`Protein Product`, ignore.case = TRUE),]
rifin <- kaust_var_genes[grep("rifin",kaust_var_genes$`Protein Product`, ignore.case = TRUE),]
pfmc2TM<- kaust_var_genes[grep("Pfmc-2TM",kaust_var_genes$`Protein Product`,ignore.case = T),]
var_genes <- kaust_var_genes[grep("PfEMP1",kaust_var_genes$`Protein Product`,ignore.case = T),]
source('~/rohit/kaust_proj/sample_swap_corrected/all_functions.R')
############ Data Pre-processing Chunk ####################


################### Header Chunk #########################

header <- dashboardHeader(title = "P. Falciparum SingleCell Analysis")
################### Sidebar Chunk ########################
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Internal QC", tabName = "iqc", icon = icon("database")),
    menuItem("Seurat QC", tabName = "sqc", icon = icon("info-circle")),
    menuItem("Integration with MCA", tabName = "mcaqc", icon = icon("info-circle")),
    menuItem("Var Gene Analysis", tabName = "vgene", icon = icon("info-circle")),
    menuItem("GO Enrichment (DE Genes)", tabName = "go_enrich", icon = icon("info-circle")),
    menuItem("GO Enrichment (Conserved Markers)", tabName = "go_enrich_conserved", icon = icon("info-circle"))
    )
  )


################### Body Chunk ##########################
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "iqc", 
            selectInput("smp","Select a Sample ",c("16D","16R","40D","40R")),
            imageOutput("pt")
            ),
    tabItem(tabName = "sqc",
            selectInput("gpby","Select x.axis Group For Vilon Plot",c("Sample","batch","seurat_clusters")),
            box(plotOutput("vlnplt"), solidHeader = T, collapsible = T, title = "QC Vilon Plot", status = "primary"),
            box(plotOutput("featureplt"), solidHeader = T, collapsible = T, title = "QC Featue-Feature Plot", status = "primary")
    ),
    tabItem(tabName = "mcaqc" ,
            box(plotOutput("p3"), solidHeader = T, width = 8, collapsible = T, title = "Samplewise", status = "primary"),
            box(plotOutput("p1"), solidHeader = T, width = 4, collapsible = T, title = "Overlap with MCA Data", status = "primary"),
            box(plotOutput("p2"), solidHeader = T, width = 4, collapsible = T, title = "Batchwise", status = "primary"),
            box(plotOutput("p4"), solidHeader = T, width = 4, collapsible = T, title = "Clusterwise", status = "primary")
            ),
    
    tabItem(tabName = "vgene",
            selectInput("vg1","Select Var-Gene",unique(var_genes$`Gene ID`)),
            box(plotOutput("vg2"), solidHeader = T,width = 12, collapsible = T, title = "Var-Gene Expression", status = "primary"),
            selectInput("vg3","Select Per Cell Expression",
                        list("Var Gene Exp"="var_genes_percent","Stevor Gene Exp" ="stevor_percent","Rifin Gene Exp"="rifin_percent","Surfin Gene Exp"="surfin_percent","PFC2TM Gene Exp"="pfc2tm_percent")),
            box(plotOutput("vg4"), solidHeader = T, width = 12,collapsible = T, title = "Gene Expression", status = "primary"),
            box(plotOutput("vg5"), solidHeader = T,collapsible = T, title = "Gene Expression", status = "primary"),
            box(plotOutput("vg6"), solidHeader = T,collapsible = T, title = "Number of Var Gene Per Cell", status = "primary")
            ),
    tabItem(tabName = "go_enrich",
            selectInput("mrkgenes_16hpi","Select the Cluster for Marker Genes (16D vs 16R)",paste("Cluster",clusters_16hpi)),
            selectInput("mrkgenes_40hpi","Select the Cluster for Marker Genes (40D vs 40R)",paste("Cluster",clusters_40hpi)),      box(plotOutput("dotplt_16hpi"), solidHeader = T,width=5, collapsible = T, title = "Gene Ontology (16D vs 16R)", status = "primary"),
            box(plotOutput("dotplt_40hpi"), solidHeader = T,width=5, collapsible = T, title = "Gene Ontology (40D vs 40R)", status = "primary"),
            box(plotOutput("barplt_16hpi"), solidHeader = T,width=5,collapsible = T, title = "GO Enriched Biological Pathways (16D vs 16R)", status = "primary"),
            box(plotOutput("barplt_40hpi"), solidHeader = T,width=5,collapsible = T, title = "GO Enriched Biological Pathways (40D vs 40R)", status = "primary"),
            box(plotOutput("wordplt_16hpi"), solidHeader = T,width=5,collapsible = T, title = "WordCloud (16D vs 16R)", status = "primary"),
            box(plotOutput("wordplt_40hpi"), solidHeader = T,width=5,collapsible = T, title = "WordCloud (40D vs 40R)", status = "primary")
    ),
    tabItem(tabName = "go_enrich_conserved",
            selectInput("mrkgenes_cons","Select the Cluster for Marker Genes (40D vs 40R)",paste("Cluster",clusters_40hpi)),
            box(plotOutput("dotplt_cons"), solidHeader = T,width=5, collapsible = T, title = "Gene Ontology (conserved markers)", status = "primary"),
            box(plotOutput("barplt_cons"), solidHeader = T,width=5,collapsible = T, title = "GO Enriched Biological Pathways (conserved markers)", status = "primary"),
            box(plotOutput("wordplt_cons"), solidHeader = T,width=5,collapsible = T, title = "WordCloud (conserved markers)", status = "primary")
  )))


############ ui chunk ####################################
ui <- dashboardPage(skin = "black", header, sidebar, body)

############ Server Chunk ################################
server <- function(input, output, session,...) {
  cdata <- session$clientData
##---------------------------Tab1--------------------------##
    output$pt <- renderImage({filename <- normalizePath(file.path('/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/www',
                                        paste0('bcrank_p', input$smp, '.png')))# Return a list containining the filename
  list(src = filename)}, deleteFile = FALSE)
 

##--------------------------Tab2----------------------------##

temp <- seurat.obj[,seurat.obj@meta.data$Sample=="16D"|seurat.obj@meta.data$Sample=="16R"|seurat.obj@meta.data$Sample=="40D"|seurat.obj@meta.data$Sample=="40R"]
output$vlnplt <- renderPlot(VlnPlot(temp,pt.size=-1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = input$gpby))
mydata <- reactive({s <- input$gpby})

output$featureplt <- renderPlot({
  s <- mydata()
  p1 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = s)
  p2 <- FeatureScatter(temp, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = s)
  p1+p2})

##----------------------Tab 3------------------------------##

output$p1 <- renderPlot(DimPlot(seurat.obj, reduction = "umap", label = TRUE, label.size = 8, repel = TRUE,group.by = "bulk",cols=c('early_troph' = '#72147e', 'late_troph' = 'green', 'ring' = '#04009a', 'schizont'='#fa9905', 'Kaust Single Cells'='grey'))+size)
output$p2 <- renderPlot(DimPlot(seurat.obj, reduction = "umap", group.by = "batch",label.size = 8,cols = c('16hpi' = '#72147e', '40hpi' = '#fa9905', 'MCA'='grey'))+size)
output$p4 <- renderPlot(DimPlot(seurat.obj, reduction = "umap", label = TRUE, label.size = 8,repel = TRUE)+size)
output$p3 <- renderPlot(DimPlot(seurat.obj, reduction = "umap", group.by = "Sample",label.size = 8, split.by = "batch",cols = c('16D' = '#72147e', '16R' = 'green', '40D' = '#04009a', '40R' = '#ff5200', 'MCA'='grey'))+size)

##-------------------Tab 4 ---------------------------------##
output$vg2 <- renderPlot(FeaturePlot(seurat.obj,features = input$vg1,split.by = "Sample",by.col = T,label = T,label.size = 6))
s <- reactive({features <- input$vg3})
output$vg4 <- renderPlot(FeaturePlot(seurat.obj,features = input$vg3,split.by = "Sample",by.col = T,label = T,label.size = 6))
output$vg5 <- renderPlot({
features <- s()
df <- data.frame(seurat.obj@meta.data[,features],seurat.obj@meta.data$Sample,seurat.obj$batch)
colnames(df) <- c("perc","sample","batch")
df$group <- df$batch
df <- df[df$batch==c("16hpi","40hpi"),]
for(i in 1:nrow(df)){
  ifelse(df$sample[i]=="16D"| df$sample[i]=="40D",df$group[i] <- "Normal",df$group[i] <- "Knockout") 
}
ggboxplot(df, x="batch",y="perc", color = "group", outlier.shape = NA,palette = "jco")+stat_compare_means(aes(group=group))+xlab("Batch")+ylab("Average Gene Expression")+size})
df <- data.frame(table(seurat.obj$Sample,seurat.obj$var_num))
output$vg6 <- renderPlot(ggplot(df, aes(fill=Var1, y=Freq, x=Var2)) + 
                           geom_bar(position="dodge", stat="identity")+theme_bw()+scale_fill_viridis(discrete = T)+xlab("No. of Var Genes")+ylab("No. of Cells")+labs(fill="Sample")+size)

##---------------Tab5----------------------------------------##
list.of.files_16hpi <- as.list(paste0("de_",clusters_16hpi,"_16hpi.csv"))
names(list.of.files_16hpi) <- paste("Cluster",clusters_16hpi)
file_16hpi <- reactive({cluster <- list.of.files_16hpi[[input$mrkgenes_16hpi]]})

list.of.files_40hpi <- as.list(paste0("de_",clusters_40hpi,"_40hpi.csv"))
names(list.of.files_40hpi) <- paste("Cluster",clusters_40hpi)
file_40hpi <- reactive({cluster <- list.of.files_40hpi[[input$mrkgenes_40hpi]]})

output$dotplt_16hpi <- renderPlot({
  cluster <- file_16hpi()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
go <- enrichGO(gene=genes, 
               ont ="ALL",keyType="SYMBOL",
               OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
               qvalueCutoff = 0.10, 
               pAdjustMethod = "BH")
dotplot(go, showCategory=10)+size
})

output$barplt_16hpi <- renderPlot({
  cluster <- file_16hpi()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
                 qvalueCutoff = 0.10, 
                 pAdjustMethod = "BH")
  barplot(go, 
        drop = TRUE, 
        showCategory = 10,
        font.size = 8)+size})
#View(as.data.frame(go))
## Frequency
## Some clusters were removed because they had no cluster specific marker
## enriched in GO Term assessment (For example Cluster 10)
output$wordplt_16hpi <- renderPlot({
  cluster <- file_16hpi()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.10, minGSSize = 3,
                 pAdjustMethod = "BH")
  dotplot(go, showCategory=10)
  wcdf<-read.table(text=go$GeneRatio, sep = "/")[1]
## function terms
wcdf$term<-go[,3]
par(mar = rep(0, 4))
wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
})

output$dotplt_40hpi <- renderPlot({
  cluster <- file_40hpi()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
                 qvalueCutoff = 0.10, 
                 pAdjustMethod = "BH")
  dotplot(go, showCategory=10)+size
})

output$barplt_40hpi <- renderPlot({
  cluster <- file_40hpi()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
                 qvalueCutoff = 0.10, 
                 pAdjustMethod = "BH")
  barplot(go, 
          drop = TRUE, 
          showCategory = 10,
          font.size = 8)+size})
#View(as.data.frame(go))
## Frequency
## Some clusters were removed because they had no cluster specific marker
## enriched in GO Term assessment (For example Cluster 10)
output$wordplt_40hpi <- renderPlot({
  cluster <- file_40hpi()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.10, minGSSize = 3,
                 pAdjustMethod = "BH")
  dotplot(go, showCategory=10)
  wcdf<-read.table(text=go$GeneRatio, sep = "/")[1]
  ## function terms
  wcdf$term<-go[,3]
  par(mar = rep(0, 4))
  wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
})
###---------------------------Tab 6----------------------------##
list.of.files_cons <- as.list(paste0("cluster.",clusters_16hpi,".conserved.markers.csv"))
names(list.of.files_cons) <- paste("Cluster",clusters_16hpi)
file_cons <- reactive({cluster <- list.of.files_cons[[input$mrkgenes_cons]]})

output$dotplt_cons <- renderPlot({
  cluster <- file_cons()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
                 qvalueCutoff = 0.10, 
                 pAdjustMethod = "BH")
  dotplot(go, showCategory=10)+size
})

output$barplt_cons <- renderPlot({
  cluster <- file_cons()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, minGSSize = 3,
                 qvalueCutoff = 0.10, 
                 pAdjustMethod = "BH")
  barplot(go, 
          drop = TRUE, 
          showCategory = 10,
          font.size = 8)+size})
#View(as.data.frame(go))
## Frequency
## Some clusters were removed because they had no cluster specific marker
## enriched in GO Term assessment (For example Cluster 10)
output$wordplt_cons <- renderPlot({
  cluster <- file_cons()
  cluster.file <- read.csv(paste0("/home/pichkari/rohit/kaust_proj/sample_swap_corrected/app/",cluster))
  genes <- gsub("-","_",cluster.file$X)
  go <- enrichGO(gene=genes, 
                 ont ="ALL",keyType="SYMBOL",
                 OrgDb = pf.go,pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.10, minGSSize = 3,
                 pAdjustMethod = "BH")
  dotplot(go, showCategory=10)
  wcdf<-read.table(text=go$GeneRatio, sep = "/")[1]
  ## function terms
  wcdf$term<-go[,3]
  par(mar = rep(0, 4))
  wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(4, .1)), colors=brewer.pal(8, "Dark2"), max.words = 25)
})


}


shinyApp(ui, server)
