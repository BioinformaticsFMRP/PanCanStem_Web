library(shinydashboard)
library(shiny)
library(shinyjs)
library(shinyBS)
library(plotly)
pdf(NULL)

getTCGAdisease <- function(){
  projects <- TCGAbiolinks:::getGDCprojects()
  disease <-  projects$project_id
  idx <- grep("disease_type",colnames(projects))
  names(disease) <-  paste0(projects[[idx]], " (",disease,")")
  disease <- disease[sort(names(disease))]
  return(disease)
}

#' busyIndicator
#'
# This is a function to indicate the work is in progress, it was created for the plots
# that rendering were taking long and withprogress was not working.
# @param text The text to show
# @param wait The amount of time to wait before showing the busy indicator. The
#   default is 1000 which is 1 second.
#
# @export
busyIndicator <- function(text = "Working in progress...") {
  div(
    id = 'busyModal', class = 'modal', role = 'dialog', 'data-backdrop' = 'static',
    div(
      class = 'modal-dialog modal-sm',
      div(id = 'modal-content-busy',
          class = 'modal-content',
          div(class = 'modal-header', h4(class = 'modal-title', text)),
          div(class = 'modal-body', p(h2(HTML('<i class="fa fa-cog fa-spin"></i>'))))
      )
    )
  )
}


header <- dashboardHeader(
  title = " PanCanStem Project: Enrichment analysis for DNAss and RNAss acrosss Clinical and Molecular phenotypes",
  titleWidth = 1000
)

body <- dashboardBody(
  tagList(
    singleton(tags$head(tags$script(
      singleton(tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "pancan.css"))),
      singleton(tags$head(singleton(tags$script(src = 'events.js')))),
      singleton(tags$head(tags$style(HTML('
        .skin-blue .main-header .logo {
          background-color: #3c8dbc;
        }
        .skin-blue .main-header .logo:hover {
          background-color: #3c8dbc;
        }
      '))))
    )))),
  fluidRow(
    column(width = 9,
           bsAlert("message"),
           box(width = NULL, solidHeader = TRUE,
               DT::dataTableOutput('tbl')
           ),
           box(width = NULL, solidHeader = TRUE,
               title = "Enrichment analysis for DNAss and RNAss acrosss Clinical and Molecular phenotypes (Mutation/Molecular subtypes)",
               plotlyOutput("butterflyPlot")
           ),
           box(width = NULL, solidHeader = TRUE,
               title = "DNA",
               plotOutput("plotGseaTableDNA"),
               plotOutput("plotEnrichmentDNA")
           ),
           box(width = NULL, solidHeader = TRUE,
               title = "RNA",
               plotOutput("plotGseaTableRNA"),
               plotOutput("plotEnrichmentRNA")
           )
    ),
    column(width = 3,
           box(width = NULL, 
               selectizeInput('cancertype',
                              'Cancer type',
                              getTCGAdisease(),
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
      body),
    busyIndicator() # Add rendering in progress...
  )
)
