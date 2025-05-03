library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("VCFComparisonModule", "vcfSampleComparison.R"))

# Serverová časť modulu pre VCF porovnanie
vcfComparisonServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(TUTUUTU = "~/Documents/BP/data", Home = path.expand("~"))
    shinyDirChoose(input, "vcf_comparison_folder", roots = volumes, session = session)
    
    vcf_comparison_folder_path <- reactive({
      req(typeof(input$vcf_comparison_folder) != "integer")
      parseDirPath(volumes, input$vcf_comparison_folder)
    })
    output$selected_comparison_path <- renderText({
      folder_path <- vcf_comparison_folder_path()
      req(folder_path)
      paste("Selected folder: ",folder_path)
    })
    processed_data <- reactive({
      folder <- vcf_comparison_folder_path()
      req(folder)
      withProgress(message = "Preprocessing files...", value = 0, {
        prepare_vcf_files(folder)
      })
    })
    filtered_data <- reactive({
      data <- processed_data()
      req(data)
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      data
    })
    observe({
      data <- processed_data()
      req(data)
      families <- c("All", unique(data$family))
      updateSelectInput(session, "family_select", choices = as.character(families), selected = families[1])
    })
    observe({
      data <- processed_data()
      req(data)
      files <- c("-" ,unique(data$name))
      updateSelectInput(session, "first_file_select", choices = as.character(files))
    })
    observe({
      data <- processed_data()
      req(data)
      files <- c("-" ,unique(data$name))
      updateSelectInput(session, "second_file_select", choices = as.character(files))
    })
    
    # Mutation Counts
    output$num_of_mutation <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      num_of_mutation(data)
    })
    output$mutation_heatmap <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      mutation_heatmap(data)
    })
    output$mutation_types_distribution <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      mutation_types_distribution(data)
    })
    output$mutation_subtypes_distribution <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      mutation_subtypes_distribution(data)
    })
    
    # SNP Analysis
    output$snp_class_comparison <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      snp_class_comparison(data)
    })
    output$transversion_transitions <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      transversion_transitions(data)
    })
    
    # INDEL Analysis
    output$insertion_deletions <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      insertion_deletions(data)
    })
    output$indel_len_boxplot <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      indel_len_boxplot(data)
    })
    
    # Quality Analysis
    output$quality_boxplot <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      quality_boxplot(data)
    })
    
    # Allele Frequency Analysis
    output$frequency_ridges <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      frequency_ridges(data)
    })
    
    # Read Depth Analysis
    output$read_depth <- renderPlot({
      data <- processed_data()
      if (input$family_select != "All") {data %<>% filter(family == input$family_select)}
      read_depth(data)
    })
    
    # Venn diagram
    output$venn_diagram <- renderPlot({
      data <- processed_data()
      file1 <- input$first_file_select
      file2 <- input$second_file_select
      req(file1 != "-")
      req(file2 != "-")
      venn_diagram(file1, file2, data)
    })
    
    
    # PDF REPORT
    output$download_vcf_comparison_pdf <- downloadHandler(
      filename ="vcfComparisonReport.pdf",
      content = function(file) {
        withProgress(message = "Generating VCF comparison report...", value = 0, {
          temp_report <- file.path("VCFComparisonModule", "vcfComparisonReport.Rmd")
          file.copy("vcfComparisonReport.Rmd", temp_report, overwrite = TRUE)
          incProgress(0.2, detail = "Loading comparison data...")
          params <- list(
            folder_name = vcf_comparison_folder_path(),
            selected_family = input$family_select,
            num_of_mutation = num_of_mutation(processed_data()),
            mutation_heatmap = mutation_heatmap(processed_data()),
            mutation_types_distribution = mutation_types_distribution(processed_data()),
            mutation_subtypes_distribution = mutation_subtypes_distribution(processed_data()),
            snp_class_comparison = snp_class_comparison(processed_data()),
            transversion_transitions = transversion_transitions(processed_data()),
            insertion_deletions = insertion_deletions(processed_data()),
            indel_len_boxplot = indel_len_boxplot(processed_data()),
            quality_boxplot = quality_boxplot(processed_data()),
            frequency_ridges = frequency_ridges(processed_data()),
            read_depth = read_depth(processed_data())
          )
          incProgress(0.5, detail = "Rendering report...")
          rmarkdown::render(
            temp_report,
            output_file = file,
            params = params,
            envir = new.env(parent = globalenv())
          )
          incProgress(1, detail = "Report completed!")
        })
      }
    )
    
    
    # COMBINED VCF COMPARISON
    output$vcf_module_compare <- renderUI({
      req(processed_data())
      if (length(processed_data())==0){
        stop("Folder does not contain any VCF files.")
      }
      ns <- NS(id)
      tagList(
        tags$h5(
          "Filters",
          tags$span(
            icon("info-circle"),
            title = "Focuses the comparative analysis on samples from a selected family, allowing insight into mutation patterns across related individuals",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, selectInput(inputId = ns("family_select"), 
                                label = "Select Family:", 
                                choices = c(""), 
                                selected = ""))
        ),
        tags$h4(
          "Mutation Analysis",
          tags$span(
            icon("info-circle"),
            title = "Summarizes mutation types and their distribution to detect mutation hotspots",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, plotOutput(ns("num_of_mutation"))),
          column(6, plotOutput(ns("mutation_heatmap"))),
          tags$hr(),
          column(6, plotOutput(ns("mutation_types_distribution"))),
          column(6, plotOutput(ns("mutation_subtypes_distribution")))
        ),
        tags$hr(),
        tags$h4(
          "Single Nucleotide Polymorphism (SNP) Analysis",
          tags$span(
            icon("info-circle"),
            title = "Shows detailed view of SNP types and classifications",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
            column(6, plotOutput(ns("snp_class_comparison"))),
            column(6, plotOutput(ns("transversion_transitions")))
        ),
        tags$hr(),
        tags$h4(
          "Insertion and Deletion (INDEL) Analysis",
          tags$span(
            icon("info-circle"),
            title = "Shows analysis of insertions and deletions, including length statistics",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
            column(6, plotOutput(ns("insertion_deletions"))),
            column(6, plotOutput(ns("indel_len_boxplot")))
        ),
        tags$hr(),
        tags$h4(
          "Quality Analysis",
          tags$span(
            icon("info-circle"),
            title = "Displays the distribution of quality scores across files",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(12, plotOutput(ns("quality_boxplot")))
        ),
        tags$hr(),
        tags$h4(
          "Allele Frequency Analysis",
          tags$span(
            icon("info-circle"),
            title = "Displays the distribution of allele frequencies across files",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(12, plotOutput(ns("frequency_ridges")))
        ),
        tags$hr(),
        tags$h4(
          "Read Depth Analysis",
          tags$span(
            icon("info-circle"),
            title = "Displays the distribution of read depth across files",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(12, plotOutput(ns("read_depth")))
        ),
        tags$h4(
          "Venn Diagrams",
          tags$span(
            icon("info-circle"),
            title = "",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        tags$h5("Choose Files"),
        fluidRow(
          column(6, selectInput(inputId = ns("first_file_select"), 
                                label = "Select First File:", 
                                choices = c(""), 
                                selected = ""))
        ),
        fluidRow(
          column(6, selectInput(inputId = ns("second_file_select"), 
                                label = "Select Second File:", 
                                choices = c(""), 
                                selected = ""))
        ),
        fluidRow(column(12, plotOutput(ns("venn_diagram")))),
        downloadButton(ns("download_vcf_comparison_pdf"), "Download PDF Report"),
      )
    })
    
  })
}

