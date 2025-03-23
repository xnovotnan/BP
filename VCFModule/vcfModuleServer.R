library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)
source("data_processing.R")
source(file.path("VCFModule", "mutationAnalysis.R"))
options(shiny.maxRequestSize = 2000 * 1024^2)


# Serverová časť modulu pre Qualimap analýzu
vcfModuleServer <- function(id, processedData, chromSelectVal) {
  moduleServer(id, function(input, output, session) {
    
    # Mutation counts
    output$mut_summary <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if (chromSelect != "All") {data %<>% filter(CHROM == chromSelect)}
      
      mut_summary(data, input$analysisLevel == "Subtypes", input$valueType)
    })
    output$mut_dist <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      mut_dist(data, input$analysisLevel == "Subtypes", input$valueType)
    })
    
    # SNV Analysis
    output$snv_types <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      SNV_types(data)
    })
    output$snv_class <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      SNV_class(data)
    })
    output$snv_class_combined <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snv_class_combined(data)
    })
    output$snv_class_stacked <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snv_class_stacked(data)
    })

  })
}
