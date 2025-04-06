library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("VCFComparisonModule", "vcfSampleComparison.R"))

# Serverová časť modulu pre VCF porovnanie
vcfComparisonServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(TUTUUTU = "~/Documents/BP/data", Home = path.expand("~"))
    shinyDirChoose(input, "vcfComparisonFolder", roots = volumes, session = session)
    
    vcfComparisonFolderPath <- reactive({
      req(typeof(input$vcfComparisonFolder) != "integer")
      parseDirPath(volumes, input$vcfComparisonFolder)
    })
    output$selectedComparisonPath <- renderText({
      folderPath <- vcfComparisonFolderPath()
      req(folderPath)
      paste("Selected folder: ",folderPath)
    })
    processedData <- reactive({
      folder <- vcfComparisonFolderPath()
      req(folder)
      prepare_vcf_files(folder)
    })
    observe({
      data <- processedData()
      req(data)
      families <- c("All", unique(data$family))
      updateSelectInput(session, "family_select", choices = as.character(families), selected = families[1])
    })
    
    # Mutation Counts
    output$num_of_mutation <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      num_of_mutation(data)
    })
    output$mutation_heatmap <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      mutation_heatmap(data)
    })
    output$mutation_types_distribution <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      mutation_types_distribution(data)
    })
    output$mutation_subtypes_distribution <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      mutation_subtypes_distribution(data)
    })
    
    # SNP Analysis
    output$snp_class_comparison <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      snp_class_comparison(data)
    })
    output$transversion_transitions <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      transversion_transitions(data)
    })
    
    # INDEL Analysis
    output$insertion_deletions <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      insertion_deletions(data)
    })
    output$indel_len_boxplot <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      indel_len_boxplot(data)
    })
    
    # Quality Analysis
    output$quality_boxplot <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      quality_boxplot(data)
    })
    
    # Allele Frequency Analysis
    output$frequency_ridges <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      frequency_ridges(data)
    })
    
    # Read Depth Analysis
    output$read_depth <- renderPlot({
      data <- processedData()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      read_depth(data)
    })
    
    # PDF REPORT
    output$download_vcf_comparison_pdf <- downloadHandler(
      filename = "vcfComparisonReport.pdf",
      content = function(file) {
        tempReport <- file.path("VCFComparisonModule", "vcfComparisonReport.Rmd")
        file.copy("vcfComparisonReport.Rmd", tempReport, overwrite = TRUE)
        params <- list(
          folder_name = input$vcfComparisonFolder,
          selected_family = input$family_select,
          num_of_mutation = num_of_mutation(processedData()),
          mutation_heatmap = mutation_heatmap(processedData()),
          mutation_types_distribution = mutation_types_distribution(processedData()),
          mutation_subtypes_distribution = mutation_subtypes_distribution(processedData()),
          snp_class_comparison = snp_class_comparison(processedData()),
          transversion_transitions = transversion_transitions(processedData()),
          insertion_deletions = insertion_deletions(processedData()),
          indel_len_boxplot = indel_len_boxplot(processedData()),
          quality_boxplot = quality_boxplot(processedData()),
          frequency_ridges = frequency_ridges(processedData()),
          read_depth = read_depth(processedData())
        )
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv()))
      }
    )
    
    
    # COMBINED VCF COMPARISON
    output$vcfModuleCompare <- renderUI({
      req(processedData())
      ns <- NS(id)
      tagList(
        h5("Filters"),
        fluidRow(
          column(6, selectInput(inputId = ns("family_select"), 
                                label = "Select Family:", 
                                choices = c(""), 
                                selected = ""))
        ),
        h4("Mutation counts"),
        fluidRow(
          column(6, plotOutput(ns("num_of_mutation"))),
          column(6, plotOutput(ns("mutation_heatmap"))),
          tags$hr(),
          column(6, plotOutput(ns("mutation_types_distribution"))),
          column(6, plotOutput(ns("mutation_subtypes_distribution")))
        ),
        tags$hr(),
        h4("Single Nucleotide Polymorphism (SNP) Analysis"),
        fluidRow(
            column(6, plotOutput(ns("snp_class_comparison"))),
            column(6, plotOutput(ns("transversion_transitions")))
        ),
        tags$hr(),
        h4("Insertion and Deletion (INDEL) Analysis"),
        fluidRow(
            column(6, plotOutput(ns("insertion_deletions"))),
            column(6, plotOutput(ns("indel_len_boxplot")))
        ),
        tags$hr(),
        h4("Quality Analysis"),
        fluidRow(
          column(12, plotOutput(ns("quality_boxplot")))
        ),
        tags$hr(),
        h4("Allele Frequency Analysis"),
        fluidRow(
          column(12, plotOutput(ns("frequency_ridges")))
        ),
        tags$hr(),
        h4("Read Depth Analysis"),
        fluidRow(
          column(12, plotOutput(ns("read_depth")))
        ),
        downloadButton(ns("download_vcf_comparison_pdf"), "Download PDF Report")
      )
    })
  })
}

