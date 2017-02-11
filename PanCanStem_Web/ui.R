library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinyBS)
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

header <- dashboardHeader(
  title = "PanCanStem_Web"
)

body <- dashboardBody(
  tagList(
    singleton(tags$head(tags$script(
      singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "pancan.css")))
    )))),
  fluidRow(
    column(width = 9,
           bsAlert("message"),
           box(width = NULL, solidHeader = TRUE,
               plotOutput("plotGseaTable")
           ),
           box(width = NULL,
               plotOutput("distPlot")
           )
    ),
    column(width = 3,
           box(width = NULL, 
               selectizeInput('experiment',
                              'Experiment filter',
                              c("Gene expression","DNA methylation"),
                              multiple = FALSE),
               selectizeInput('cancertype',
                              'Cancer type',
                              NULL,
                              #getTCGAdisease(),
                              multiple = FALSE),
               selectizeInput('feature',
                              'Feature',
                              NULL,
                              multiple = FALSE),
               actionButton("calculate",
                            "Calculate",
                            style = "background-color: #000080;
                            color: #FFFFFF;
                            margin-left: auto;
                            margin-right: auto;
                            width: 100%",
                            icon = icon("flask"))
           ),
           box(width = NULL, 
               selectizeInput('pathway',
                              'Pathway',
                              NULL,
                              multiple = FALSE),
               actionButton("plot",
                            "Plot GSEA enrichment",
                            style = "background-color: #000080;
                            color: #FFFFFF;
                            margin-left: auto;
                            margin-right: auto;
                            width: 100%",
                            icon = icon("flask"))
           )
    )
  )
)

shinyUI(
  bootstrapPage(
    useShinyjs(),
    div(id = "loading-content",
        img(src = "loading.gif")
        ),
    dashboardPage(
      skin = "blue",
      header,
      dashboardSidebar(disable = TRUE),
      body)
  )
)
