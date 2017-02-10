#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

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
      if(input$experiment == "Gene expression") ret <- as.character(sort(unique(pd.maf.450$cancer.type)))
      if(input$experiment == "DNA methylation") ret <- as.character(sort(unique(pd.maf.RNA$cancer.type)))
      ret
    }, server = TRUE)
  })
  
  volcano.values <- reactive({
    feature <- "IDH1"
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
    test <- subset(pd, cancer.type %in% levels(pd$cancer.type) & TCGAlong.id %in% primary) 
    test <- droplevels(test)
    stats <- test[order(test[,col]),] #rank samples
    stats <- structure(stats[,col], names=rownames(stats))
    pathways <- as.list(unstack(test[,c("TCGAlong.id", feature )]))
    withProgress(message = 'Creating plot',
                 detail = 'This may take a while...', value = 0, {
    result <- fgsea(pathways = pathways, stats = stats,  nperm=10, minSize=5, maxSize=500)
    })
    ret <- list(stats = stats, 
                pathways = pathways)
    return(result)
  })
  
  observeEvent(input$plot , {
    print("Hey2")
    output$distPlot <- renderPlot({
     print("Hey")
     ret <- volcano.values()
     plotEnrichment(ret$pathways[["Missense_Mutation"]], ret$stats) + labs(title="IDH1-Mutant")
     #plotEnrichment(pathways[["WT"]], stats) + labs(title="IDH1-WT")
     #plotGseaTable(pathways, stats, result,  gseaParam = 0.5)
   })
  })
  
})
