#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

getTCGAdisease <- function(){
  projects <- TCGAbiolinks:::getGDCprojects()
  disease <-  projects$project_id
  names(disease) <-  paste0(projects$disease_type, " (",disease,")")
  disease <- disease[sort(names(disease))]
  tcga.disease <- disease[grep("TCGA",disease)]
  tcga.disease <- gsub("TCGA-","",tcga.disease)
  tcga.disease <- tcga.disease[sort(names(tcga.disease), index.return=TRUE)$ix]
  return(tcga.disease)
}
# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("PanCanStem_Web"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      selectizeInput('experiment',
                     'Experiment filter',
                     c("Gene expression","DNA methylation"),
                     multiple = FALSE),
      selectizeInput('cancertype',
                     'Cancer type',
                     getTCGAdisease(),
                     multiple = FALSE),
      selectizeInput('feature',
                     'Feature',
                     NULL,
                     multiple = FALSE),
      selectizeInput('featureLevels',
                     'Feature Levels',
                     NULL,
                     multiple = FALSE),
      actionButton("plot",
                   "Plot",
                   style = "background-color: #000080;
                                    color: #FFFFFF;
                                    margin-left: auto;
                                    margin-right: auto;
                                    width: 100%",
                   icon = icon("flask"))
    ),
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
))
