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
    
    output$mut_summary <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if (chromSelect != "All") {data %<>% filter(CHROM == chromSelect)}
      
      mut_summary(data, input$analysisLevel == "Subtypes")
    })
    
    output$mut_dist <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      mut_dist(data, input$analysisLevel == "Subtypes", input$valueType)
    })
    
    
  })
}
