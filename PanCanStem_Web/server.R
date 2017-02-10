#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#
library(data.table)
library(shiny)
library(ggplot2)
library(fgsea)
load("data/pd.450.prim_20170207.Rda")
load("data/pd.maf.450.Rda")
load("data/pd.maf.RNA.Rda")
load("data/pd.mRNA.prim_20170208.Rda")
set.seed(10)



# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  
  
  observeEvent(input$experiment, {
    # This must be checked how do we have smoking?
    updateSelectizeInput(session, 'feature', choices = {
      if(isolate({input$experiment}) == "Gene expression") ret <- as.character(colnames(pd.maf.RNA))
      if(isolate({input$experiment}) == "DNA methylation")  ret <- as.character(colnames(pd.maf.450))
      ret
    }, server = TRUE)
  })
  
  observeEvent(input$feature, {
    # This must be checked how do we have smoking?
    updateSelectizeInput(session, 'featureLevels', choices = {
      if(is.null(input$feature) || input$feature == "") {
        ret <- NULL
      } else if(isolate({input$experiment}) == "Gene expression") {
        ret <- as.character(unique(pd.maf.RNA[,input$feature]))
      } else if(isolate({input$experiment}) == "DNA methylation") {
        ret <- as.character(unique(pd.maf.450[,input$feature]))
      }
      ret
    }, server = TRUE)
  })
  
  volcano.values <- reactive({
    feature <- input$feature
    if(isolate({input$experiment}) == "Gene expression")  {
      pd <- pd.maf.RNA
      primary <- pd.mRNA.prim[pd.mRNA.prim$sample.type %in% c("01","03"),"TCGAlong.id"]
      col <- "RNAss"
    }
    if(isolate({input$experiment})  == "DNA methylation")  {
      pd <- pd.maf.450
      primary <- pd.450.prim[pd.450.prim$sample.type %in% c("01","03"),"TCGAlong.id"]
      col <- "DNAss"
    }
    withProgress(message = 'Calculating...',
                 detail = 'This may take a while...', value = 0, {
                   test <- subset(pd, pd$cancer.type %in% cancer.type  & pd$TCGAlong.id %in% primary) 
                   test <- droplevels(test)
                   stats <- test[order(test[,col]),] #rank samples
                   stats <- structure(stats[,col], names=as.character(stats$TCGAlong.id))
                   pathways <- as.list(unstack(test[,c("TCGAlong.id", feature )]))
                   
                   result <- fgsea(pathways = pathways, stats = stats,  nperm=10, minSize=5, maxSize=500)
                 })
    ret <- list(stats = stats, 
                result = result,
                pathways = pathways)
    return(ret)
  })
  
  observeEvent(input$plot , {
    output$distPlot <- renderPlot({
      ret <- volcano.values()
      plotEnrichment(ret$pathways[[input$featureLevels]], ret$stats) + labs(title=input$featureLevels)
      # plotGseaTable(ret$pathways, ret$stats, ret$result,  gseaParam = 0.5)
    })
  })
  
})
