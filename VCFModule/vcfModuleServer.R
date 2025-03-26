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
    
    # SNP Analysis
    output$snp_values <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      values <- snp_values(data)
      result <- sprintf("SNP Count: %s (%s %%)", label_comma()(values$count), values$percentage)
      
    })
    output$snp_types <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snp_types(data)
    })
    output$snp_class <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snp_class_barplot(data)
    })
    output$snp_class_combined <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snp_class_boxplot(data)
    })
    output$snp_class_stacked <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      snp_class_stacked(data)
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
    
    # Quality
    output$qual_avg <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      avg <- quality_avg(data)
      sprintf("Mean Quality: %s", avg)
    })
    output$qual_med <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      med <- quality_med(data)
      sprintf("Median Quality: %s", med)
    })
    output$quality_bar <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      quality_bar(data)
    })
    output$quality_on_chroms <- renderPlot({
      data <- processedData()
      req(data)
      quality_on_chroms(data)
    })
    
    # Allele frequency
    output$allele_freq_avg <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      avg <- allele_freq_avg(data)
      sprintf("Mean Allele Frequency: %s", avg)
    })
    output$allele_freq_med <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      med <- allele_freq_med(data)
      sprintf("Median Allele Frequency: %s", med)
    })
    output$allele_freq_hexbin <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      allele_freq_hexbin(data)
    })
    output$allele_freq_on_chroms <- renderPlot({
      data <- processedData()
      req(data)
      allele_freq_on_chroms(data)
    })
      
    # Read Depth
    output$read_depth_avg <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      avg <- read_depth_avg(data)
      sprintf("Mean Read Depth: %s", avg)
    })
    output$read_depth_med <- renderText({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      med <- read_depth_med(data)
      sprintf("Median Read Depth: %s", med)
    })
    output$read_depth_density <- renderPlot({
      data <- processedData()
      req(data)
      chromSelect <- chromSelectVal()
      if(chromSelect != "All"){data %<>% filter(CHROM == chromSelect)}
      read_depth_density(data)
    })
    output$read_depth_on_chroms <- renderPlot({
      data <- processedData()
      req(data)
      read_depth_on_chroms(data)
    })
    
  })
}
