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
library(shinyBS)
load("data/pd.450.prim_20170207.Rda")
load("data/pd.maf.450.Rda")
load("data/pd.maf.RNA.Rda")
load("data/pd.mRNA.prim_20170208.Rda")
load("data/features.Rda")
set.seed(10)

# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  
  observe({
    updateSelectizeInput(session, 'feature', choices = {
      sort(unique(as.character(features)))
    }, server = TRUE)
    updateSelectizeInput(session, 'cancertype', choices = {
      sort(unique(pd.maf.450$cancer.type))
    }, server = TRUE)
  })
  
  
  observeEvent(input$feature, {
    updateSelectizeInput(session, 'pathway', choices = {
      if(is.null(input$feature) || input$feature == "") {
        ret <- NULL
      } else if(isolate({input$experiment}) == "Gene expression") {
        if(input$feature %in% colnames(pd.maf.RNA)) {
          ret <- levels(pd.maf.RNA[,get(input$feature)])
        } else {
          ret <- as.character(unique(pd.mRNA.prim[,get(input$feature)]))
        }
      } else if(isolate({input$experiment}) == "DNA methylation") {
        if(input$feature %in% colnames(pd.maf.450)) {
          ret <- levels(pd.maf.450[,get(input$feature)])
        } else {
          ret <- as.character(unique(pd.450.prim[,get(input$feature)]))
        }
      }
      ret
    }, server = TRUE)
  })
  
  volcano.values <- reactive({
    closeAlert(session, "Alert")
    if(input$calculate){
      feature <- isolate({input$feature})
      if(isolate({input$experiment}) == "Gene expression")  {
        if(feature %in% colnames(pd.maf.RNA)) {
          pd <- pd.maf.RNA
          nperm <- 1000
        } else {
          pd <- pd.mRNA.prim
          nperm <- 10000
        }
        primary <- as.character(pd.mRNA.prim[pd.mRNA.prim$sample.type %in% c("01","03"),get("TCGAlong.id")])
        col <- "RNAss"
      }
      if(isolate({input$experiment}) == "DNA methylation")  {
        if(feature %in% colnames(pd.maf.450)) {
          pd <- pd.maf.450
          nperm <- 1000
        } else {
          pd <- pd.450.prim
          nperm <- 10000
        }
        primary <- as.character(pd.450.prim[pd.450.prim$sample.type %in% c("01","03"),get("TCGAlong.id")])
        col <- "DNAss"
      }
      progress <- shiny::Progress$new()
      progress$set(message = "Calculating", value = 0, detail = "Preparing data")
      on.exit(progress$close())
      test <- subset(pd, pd$cancer.type %in% isolate({input$cancertype})  & pd$TCGAlong.id %in% primary) 
      stats <- test[order(test[,col,with=FALSE]),] #rank samples
      stats <- structure(stats[,get(col)], names = as.character(stats$TCGAlong.id))
      if(length(unique(test[,get(feature)])) < 2){
        createAlert(session, "message", "Alert", title = "Data input error", style =  "danger",
                    content = "There is not two levels for this combination. We cannot execute Gene Set Enrichment Analysis for this case.", append = FALSE)
        return(NULL)
      }
      pathways <- as.list(unstack(test[,c("TCGAlong.id", feature),with = FALSE]))
      progress$set(value = 0.5, detail = "Executing Gene Set Enrichment Analysis")
      result <- fgsea(pathways = pathways, stats = stats,  nperm=nperm, minSize=5, maxSize=500)
      
      ret <- list(stats = stats, 
                  result = result,
                  pathways = pathways)
      return(ret)
    }
  })
  
  observeEvent(input$plot , {
    output$distPlot <- renderPlot({
      ret <- volcano.values()
      feature.level <- isolate({input$pathway})
      if(!is.null(feature.level) & feature.level != "") 
        plotEnrichment(ret$pathways[[feature.level]], ret$stats) + labs(title=input$pathway)
    })
  })
  observeEvent(input$calculate , {
    output$plotGseaTable <- renderPlot({
      ret <- volcano.values()
      if(!is.null(ret)) plotGseaTable(ret$pathways, ret$stats, ret$result,  gseaParam = 0.5)
    })
  })
  hide("loading-content", TRUE, "fade")
})
