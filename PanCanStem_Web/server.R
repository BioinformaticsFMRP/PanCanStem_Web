library(data.table)
library(shiny)
library(ggplot2)
library(fgsea)
library(shinyBS)
library(DT)
library(plotly)

getTCGAdisease <- function(){
  projects <- TCGAbiolinks:::getGDCprojects()
  projects <- projects[grep("TCGA",projects$project_id),]
  disease <-  gsub("TCGA-","",projects$project_id)
  idx <- grep("disease_type",colnames(projects))
  names(disease) <-  paste0(projects[[idx]], " (",disease,")")
  disease <- disease[sort(names(disease))]
  return(disease)
}
load("data/pd.450.prim_20170207.Rda")
load("data/pd.maf.450.Rda")
load("data/pd.maf.RNA.Rda")
load("data/pd.mRNA.prim_20170208.Rda")
load("data/pd.all.Rda")
load("data/pd.merg.Rda")
load("data/features.Rda")
set.seed(10)

cancer.colors <- c( "ACC"="white",
                    "BLCA"="white",
                    "BRCA"="pink",
                    "CESC"="white",
                    "CHOL"="white",
                    "COAD"="white",
                    "DLBC"="white",
                    "ESCA"="white",
                    "GBM"="white",
                    "HNSC"="white",
                    "KICH"="orange",
                    "KIRC"="orange",
                    "KIRP"="orange",
                    "LAML"="white",
                    "LGG"="white",
                    "LIHC"="green",
                    "LUAD"="white",
                    "LUSC"="white",
                    "MESO"="white",
                    "OV"="	#008080", # Teal
                    "PAAD"="white",
                    "PCPG"="white",
                    "PRAD"="white",
                    "READ"="white",
                    "SARC"="white",
                    "SKCM"="white",
                    "STAD"="white",
                    "TGCT"="white",
                    "THCA"="white",
                    "THYM"="white",
                    "UCEC"="white",
                    "UCS"="white",
                    "UVM"="white"
)


# Define server logic required to draw a histogram
shinyServer(function(input, output,session) {
  
  output$tbl = DT::renderDataTable(
    datatable(pd.all,
              extensions = c('Buttons',"FixedHeader"),
              class = 'cell-border stripe',
              options = list(dom = 'Blfrtip',
                             buttons =  c('copy', 'csv', 'excel', 'pdf', 'print', 'colvis'),
                             fixedHeader = TRUE,
                             pageLength = 20,
                             scrollX = TRUE,
                             lengthMenu = list(c(10, 20, 50, -1), c('10', '20','50','All'))
              ),
              filter = 'top')  %>% 
      formatStyle(
        'padj',
        backgroundColor = styleInterval(c(0.5), c('red', 'white'))
      ) %>% 
      formatStyle(
        'cancer.type',
        backgroundColor = styleEqual(
          sort(unique(pd.all$cancer.type)), cancer.colors[sort(unique(pd.all$cancer.type))]
        )
      )
  )

  updateSelectizeInput(session, 'cancertype', choices = {
    sort(unique(as.character(getTCGAdisease())))
  }, server = TRUE)
  
  
  observeEvent(input$cancertype, {
    updateSelectizeInput(session, 'feature', choices = {
      sort(unique(as.character(pd.all[pd.all$cancer.type %in% input$cancertype,"Pt.Feature"])))
    }, server = TRUE)
  })
  
  observeEvent(input$feature, {
    feature <- isolate({input$feature})
    updateSelectizeInput(session, 'pathway', choices = {
      if(is.null(feature) || feature == "") {
        ret <- NULL
      } else if(feature %in% colnames(pd.maf.RNA)) {
        ret <- unique(levels(pd.maf.RNA[pd.maf.RNA$cancer.type %in% isolate({input$cancertype}) ,get(feature)]),
                      levels(pd.maf.450[pd.maf.450$cancer.type %in% isolate({input$cancertype}) ,get(feature)]))
      } else {
        ret <- unique(as.character(
          pd.mRNA.prim[pd.mRNA.prim$cancer.type %in% isolate({input$cancertype}),get(feature)],
          pd.450.prim[pd.450.prim$cancer.type %in% isolate({input$cancertype}),get(feature)]))
      }
      ret
    }, server = TRUE)
  })
  
  volcano.values <- reactive({
    closeAlert(session, "Alert")
    if(input$calculate){
      
      feature <- isolate({input$feature})
      if(is.null(feature) || feature == ""){ 
        createAlert(session, "message", "Alert", title = "Error", style =  "danger",
                    content = "Please select Feature", append = FALSE)
        return(NULL)
      }
      if(feature %in% colnames(pd.maf.RNA)) {
        pd <- pd.maf.RNA
        nperm <- 1000
      } else {
        pd <- pd.mRNA.prim
        nperm <- 10000
      }
      primary <- as.character(pd.mRNA.prim[pd.mRNA.prim$sample.type %in% c("01","03"),get("TCGAlong.id")])
      col <- "RNAss"
      progress <- shiny::Progress$new()
      progress$set(message = "Calculating", value = 0, detail = "Preparing data")
      on.exit(progress$close())
      test <- subset(pd, pd$cancer.type %in% isolate({input$cancertype})  & pd$TCGAlong.id %in% primary) 
      stats.rna <- test[order(test[,col,with=FALSE]),] #rank samples
      stats.rna <- structure(stats.rna[,get(col)], names = as.character(stats.rna$TCGAlong.id))
      if(length(unique(test[,get(feature)])) < 2){
        createAlert(session, "message", "Alert", title = "Data input error", style =  "danger",
                    content = "There is not two levels for this combination. We cannot execute Gene Set Enrichment Analysis for this case.", append = FALSE)
        return(NULL)
      }
      pathways.rna <- as.list(unstack(test[,c("TCGAlong.id", feature),with = FALSE]))
      progress$set(value = 0.5, detail = "Executing Gene Set Enrichment Analysis")
      result.rna <- fgsea(pathways = pathways.rna, stats = stats.rna,  nperm=nperm, minSize=5, maxSize=500)
      
      if(feature %in% colnames(pd.maf.450)) {
        pd <- pd.maf.450
      } else {
        pd <- pd.450.prim
      }
      primary <- as.character(pd.450.prim[pd.450.prim$sample.type %in% c("01","03"),get("TCGAlong.id")])
      col <- "DNAss"
      test <- subset(pd, pd$cancer.type %in% isolate({input$cancertype})  & pd$TCGAlong.id %in% primary) 
      stats.dna <- test[order(test[,col,with=FALSE]),] #rank samples
      stats.dna <- structure(stats.dna[,get(col)], names = as.character(stats.dna$TCGAlong.id))
      if(length(unique(test[,get(feature)])) < 2){
        createAlert(session, "message", "Alert", title = "Data input error", style =  "danger",
                    content = "There is not two levels for this combination. We cannot execute Gene Set Enrichment Analysis for this case.", append = FALSE)
        return(NULL)
      }
      pathways.dna <- as.list(unstack(test[,c("TCGAlong.id", feature),with = FALSE]))
      progress$set(value = 0.5, detail = "Executing Gene Set Enrichment Analysis")
      result.dna <- fgsea(pathways = pathways.dna, stats = stats.dna,  nperm=nperm, minSize=5, maxSize=500)
      
      ret <- list(stats.rna = stats.rna, 
                  stats.dna = stats.dna, 
                  result.rna = result.rna,
                  result.dna = result.dna,
                  pathways.rna = pathways.rna,
                  pathways.dna = pathways.dna)
      return(ret)
    }
  })
  
  observeEvent(input$plot , {
    output$plotEnrichmentDNA <- renderPlot({
      ret <- volcano.values()
      feature.level <- isolate({input$pathway})
      if(!is.null(feature.level) & feature.level != "") {
        if(!feature.level %in% names(ret$pathways.dna)){ 
          createAlert(session, "message", "Alert", title = "Error", style =  "danger",
                      content = "Please execute Gene Set Enrichment Analysis again.", append = FALSE)
          return(NULL)
        }
        plotEnrichment(ret$pathways.dna[[feature.level]], ret$stats.dna) + labs(title=input$pathway.dna)
      } else { 
        createAlert(session, "message", "Alert", title = "Error", style =  "danger",
                    content = "Please execute Gene Set Enrichment Analysis again.", append = FALSE)
        return(NULL)
      }
    })
  })
  
  observeEvent(input$plot , {
    output$plotEnrichmentRNA <- renderPlot({
      ret <- volcano.values()
      feature.level <- isolate({input$pathway})
      if(!is.null(feature.level) & feature.level != "") {
        if(!feature.level %in% names(ret$pathways.rna)){ 
          createAlert(session, "message", "Alert", title = "Error", style =  "danger",
                      content = "Please execute Gene Set Enrichment Analysis again.", append = FALSE)
          return(NULL)
        }
        plotEnrichment(ret$pathways.rna[[feature.level]], ret$stats.rna) + labs(title=input$pathway.rna)
      }
    })
  })
  observeEvent(input$calculate , {
    output$plotGseaTableRNA <- renderPlot({
      ret <- volcano.values()
      if(!is.null(ret)) plotGseaTable(ret$pathways.rna, ret$stats.rna, ret$result.rna,  gseaParam = 0.5,colwidths= c(1,5,1,1,1))
    })
  })
  observeEvent(input$cancertype , {
    output$butterflyPlot <- renderPlotly({
      if(!is.null(input$cancertype) & input$cancertype != "") {
        p <- ggplotly(
        ggplot(pd.merg[pd.merg$pathway.RNA %in% "Mutant" & pd.merg$cancer.type.DNA %in% input$cancertype,], 
               aes(x = NES.DNA, y = NES.RNA,tooltip = ID.RNA, color = c(padj.DNA < 0.05 | padj.RNA < 0.05))) + 
          geom_point() + 
          geom_vline(xintercept = 0) +
          geom_hline(yintercept = 0) +
          scale_color_manual(values = c("black","red")) +
          labs(colour = "Significant (padj<0.05)", 
               x = "DNAss Enrichment Score (NES)", 
               y = "RNAss Enrichment Score (NES)") +
          theme_bw() + 
          theme(legend.position="bottom")
        ) %>%   
          add_annotations( text="Significant (padj<0.05)", 
                                 xref="paper", yref="paper",
                                 x=0, xanchor="left",
                                 y=-0.2, yanchor="bottom",    # Same y as legend below
                                 legendtitle=TRUE, showarrow=FALSE ) %>%
          layout(title = paste0("DNAss vs RNAss Mutation Enrichment (",input$cancertype,")"),
                 #xaxis = list(showticklabels = FALSE),
                 legend = list(orientation = "h",
                               y = -0.2, x = 0.0))
        p$elementId <- NULL
        p
      } else {
        p <- plotly_empty(type = "scatter")
        p$elementId <- NULL
        p
      }
    })
  })
  
  
  observeEvent(input$calculate , {
    output$plotGseaTableDNA  <- renderPlot({
      ret <- volcano.values()
      if(!is.null(ret)) plotGseaTable(ret$pathways.dna, ret$stats.dna, ret$result.dna,  gseaParam = 0.5,colwidths= c(1,5,1,1,1))
    })
  })
  
  hide("loading-content", TRUE, "fade")
})


