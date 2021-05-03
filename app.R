library(shiny)
library(shinydashboard)
library(Seurat)
library(ggplot2)

############ Change Path Here #############################

path <- "C:/Users/Rohit Satyam/Desktop/cordium.csv"

############ Data Pre-processing Chunk ####################


################### Header Chunk #########################

header <- dashboardHeader(title = "KAUST PFSingleCell")
################### Sidebar Chunk ########################
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("DimPlots", tabName = "dim_plt", icon = icon("database")),
    menuItem("Feature Plot", tabName = "feature_plt", icon = icon("info-circle")),
    menuItem("Heatmap", tabName = "heat_plt", icon = icon("info-circle"))
    )
  )


################### Body Chunk ##########################
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "dim_plt",plotOutput("pt")),
    tabItem(tabName = "feature_plt", 
            selectInput("var_gene","Select a var gene ",common_kaust_var_genes),
            selectInput("split","Select a category ",c("Sample","batch")),
    plotOutput("var")),
    tabItem(tabName = "heat_plt",
            selectInput("spt","Select a category ",c("Sample","batch")),
            plotOutput("heatmap")
            )
  )
)

############ ui chunk ####################################
ui <- dashboardPage(skin = "black", header, sidebar, body)

############ Server Chunk ################################
server <- function(input, output) {
  output$pt <- renderPlot(p1+p2+p3, height = 800)
  output$var <- renderPlot(Seurat::FeaturePlot(object = plasmodium.combined, features = input$var_gene,min.cutoff = "q9", split.by = input$split  ))
  output$heatmap <- renderPlot(Seurat::DoHeatmap(plasmodium.combined, features = common_kaust_var_genes,  group.by = input$spt)+theme(axis.text.y = element_text(color="black", size=15)),height = 1100)
  }

shinyApp(ui, server)
