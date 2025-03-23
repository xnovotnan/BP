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
    output$mutation_donut <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if (chromSelect != "All") {data %<>% filter(CHROM == chromSelect)}
      
      mutation_donut(data, input$analysisLevel == "Subtypes", input$valueType)
    })
    output$mutation_distribution <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      mutation_distribution(data, input$analysisLevel == "Subtypes", input$valueType)
    })
    
    # SNV Analysis
    output$snv_values <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      values <- snv_values(data)
      result <- sprintf("SNV Count: %s (%s %%)", label_comma()(values$count), values$percentage)
      
    })
    output$snv_types <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snv_types(data)
    })
    output$snv_class <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snv_class_barplot(data)
    })
    output$snv_class_combined <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snv_class_boxplot(data)
    })
    output$snv_class_stacked <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snv_class_stacked(data)
    })
    
    # INDEL Analysis
    output$indel_values <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      values <- indel_values(data)
      result <- sprintf("INDEL Count: %s (%s %%)", label_comma()(values$count), values$percentage)
      
    })
    output$indel_types <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      indel_types(data)
    })
    output$indel_stacked <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      indel_stacked(data)
    })
    output$indel_length_avg <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      avg <- indel_length_avg(data)
      sprintf("Mean INDEL length: %s", avg)
    })
    output$indel_length_med <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      med <- indel_length_med(data)
      sprintf("Median INDEL length: %s", med)
    })
    output$indel_length <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      indel_length(data)
    })
    output$indel_length_boxplot <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      indel_length_boxplot(data)
    })
  })
}
