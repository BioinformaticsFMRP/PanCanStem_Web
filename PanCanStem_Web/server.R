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
load("data/pd.all.Rda")
load("data/pd.merge.Rda")
load("data/pd.450.prim_20170207.Rda")
load("data/pd.maf.450.Rda")
load("data/pd.maf.RNA.Rda")
load("data/pd.mRNA.prim_20170208.Rda")
set.seed(10)
cancer.type <- "ACC"
gene <- "A2ML1"
test <- subset(pd.maf.450, cancer.type %in% levels(pd.all$cancer.type) & sample.type %in% "01") #653  19705 or #672
test <- droplevels(test)
stats <- test[order(test$DNAss),] #rank samples
stats <- structure(stats$DNAss, names=rownames(stats))
pathways <- as.list(unstack(test[,c("TCGAlong.id", gene )]))
#
stats <- test[order(test$DNAss),] #rank samples
stats <- structure(stats$DNAss, names=rownames(stats))
#
result <- fgsea(pathways = pathways, stats = stats,  nperm=100, minSize=5, maxSize=500)
#
plotEnrichment(pathways[["Missense_Mutation"]], stats) + labs(title="IDH1-Mutant")
plotEnrichment(pathways[["WT"]], stats) + labs(title="IDH1-WT")
#
plotGseaTable(pathways, stats, result,  gseaParam = 0.5)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2] 
    bins <- seq(min(x), max(x), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')
    
  })
  
})
