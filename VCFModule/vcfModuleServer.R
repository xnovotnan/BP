library(shiny)
library(shinyWidgets)
library(tinytex)
source(file.path("VCFModule", "mutationAnalysis.R"))

vcfModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    processedData <- reactive({
      req(input$vcfFile)
      withProgress(message = "Preprocessing file...", value = 0, {
        prepare_data(input$vcfFile$datapath)
      })
    })
    observe({
      data <- processedData()
      req(data)
      chroms <- unique(data$CHROM)
      chromChoices <- c("All", as.character(chroms))
      updateSelectInput(session, "chrom_select", choices = chromChoices, selected = "All")
    })
    
    # Mutation counts
    output$mutation_donut <- renderPlot({
      data <- processedData()
      req(data)
      if (input$chrom_select != "All") {data %<>% filter(CHROM == input$chrom_select)}
      mutation_donut(data, input$analysisLevel == "Subtypes")
    })
    output$mutation_distribution <- renderPlot({
      data <- processedData()
      req(data)
      mutation_distribution(data, input$analysisLevel == "Subtypes")
    })
    
    # SNP Analysis
    output$snp_value <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      values <- snp_value(data)
      result <- sprintf("SNP Count: %s (%s %%)", label_comma()(values$count), values$percentage)
      
    })
    output$snp_types_donut <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_types_donut(data)
    })
    output$snp_class_stacked <- renderPlot({
      data <- processedData()
      req(data)
      snp_class_stacked(data)
    })
    output$snp_class_boxplot <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_class_boxplot(data)
    })
    output$snp_class_barplot <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_class_barplot(data)
    })
    
    # INDEL Analysis
    output$indel_values <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      values <- indel_values(data)
      result <- sprintf("INDEL Count: %s (%s %%)", label_comma()(values$count), values$percentage)
      
    })
    output$indel_types <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_types(data)
    })
    output$indel_stacked <- renderPlot({
      data <- processedData()
      req(data)
      indel_stacked(data)
    })
    output$indel_length_avg <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      avg <- indel_length_avg(data)
      sprintf("Mean INDEL length: %s", avg)
    })
    output$indel_length_med <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      med <- indel_length_med(data)
      sprintf("Median INDEL length: %s", med)
    })
    output$indel_length <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_length(data)
    })
    output$indel_length_boxplot <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_length_boxplot(data)
    })
    
    # Quality
    output$qual_avg <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      avg <- quality_avg(data)
      sprintf("Mean Quality: %s", avg)
    })
    output$qual_med <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      med <- quality_med(data)
      sprintf("Median Quality: %s", med)
    })
    output$quality_bar <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
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
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      avg <- allele_freq_avg(data)
      sprintf("Mean Allele Frequency: %s", avg)
    })
    output$allele_freq_med <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      med <- allele_freq_med(data)
      sprintf("Median Allele Frequency: %s", med)
    })
    output$allele_freq_hexbin <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
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
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      avg <- read_depth_avg(data)
      sprintf("Mean Read Depth: %s", avg)
    })
    output$read_depth_med <- renderText({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      med <- read_depth_med(data)
      sprintf("Median Read Depth: %s", med)
    })
    output$read_depth_density <- renderPlot({
      data <- processedData()
      req(data)
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      read_depth_density(data)
    })
    output$read_depth_on_chroms <- renderPlot({
      data <- processedData()
      req(data)
      read_depth_on_chroms(data)
    })
    
    
    # COMBINED VCF ANALYSIS
    output$vcfModuleAnalysis <- renderUI({
      req(processedData())
      ns <- NS(id)
      tagList(
        h5("Filters"),
        fluidRow(
          column(6, selectInput(ns("analysisLevel"), 
                                "Select Analysis Level:",
                                choices = c("Types", "Subtypes"), 
                                selected = "Types")),
          column(6, selectInput(inputId = ns("chrom_select"), 
                                label = "Select Chromosome:", 
                                choices = c("All"), 
                                selected = "All"))
        ),
        h4("Mutation counts"),
        fluidRow(
          column(6, plotOutput(ns("mutation_donut"))),
          column(6, plotOutput(ns("mutation_distribution")))
        ),
        tags$hr(),
        h4("Single Nucleotide Polymorphism (SNP) Analysis"),
        fluidRow(
          column(12, verbatimTextOutput(ns("snp_value"))),
          column(4, plotOutput(ns("snp_types_donut"))),
          column(8, plotOutput(ns("snp_class_stacked"))),
          tags$hr(),
          column(7, plotOutput(ns("snp_class_boxplot"))),
          column(5, plotOutput(ns("snp_class_barplot")))
        ),
        tags$hr(),
        h4("Insertion and Deletion (INDEL) Analysis"),
        fluidRow(
          column(12, verbatimTextOutput(ns("indel_values"))),
          column(4, plotOutput(ns("indel_types"))),
          column(8, plotOutput(ns("indel_stacked"))),
          column(6, verbatimTextOutput(ns("indel_length_avg"))),
          column(6, verbatimTextOutput(ns("indel_length_med"))),
          column(8, plotOutput(ns("indel_length"))),
          column(4, plotOutput(ns("indel_length_boxplot")))
        ),
        tags$hr(),
        h4("Quality Analysis"),
        fluidRow(
          column(6, verbatimTextOutput(ns("qual_avg"))),
          column(6, verbatimTextOutput(ns("qual_med"))),
          column(6, plotOutput(ns("quality_bar"))),
          column(6, plotOutput(ns("quality_on_chroms")))
        ),
        tags$hr(),
        h4("Allele Frequency Analysis"),
        fluidRow(
          column(6, verbatimTextOutput(ns("allele_freq_avg"))),
          column(6, verbatimTextOutput(ns("allele_freq_med"))),
          column(6, plotOutput(ns("allele_freq_hexbin"))),
          column(6, plotOutput(ns("allele_freq_on_chroms")))
        ),
        tags$hr(),
        h4("Read Depth Analysis"),
        fluidRow(
          column(6, verbatimTextOutput(ns("read_depth_avg"))),
          column(6, verbatimTextOutput(ns("read_depth_med"))),
          column(6, plotOutput(ns("read_depth_density"))),
          column(6, plotOutput(ns("read_depth_on_chroms")))
        )
      )
    })
  })
}
