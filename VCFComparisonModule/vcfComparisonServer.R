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
      updateSelectInput(session, "file1_select", choices = as.character(files))
    })
    observe({
      data <- processed_data()
      req(data)
      files <- c("-" ,unique(data$name))
      updateSelectInput(session, "file2_select", choices = as.character(files))
    })
    observe({
      data <- processed_data()
      req(data)
      files <- c("-" ,unique(data$name))
      updateSelectInput(session, "file3_select", choices = as.character(files))
    })
    
    # Mutation Counts
    num_of_mutation_plot <- reactive({
      num_of_mutation(filtered_data())
    })
    output$num_of_mutation <- renderPlot({
      num_of_mutation_plot()
    })
    mutation_heatmap_plot <- reactive({
      mutation_heatmap(filtered_data())
    })
    output$mutation_heatmap <- renderPlot({
      mutation_heatmap_plot()
    })
    mutation_types_distribution_plot <- reactive({
      mutation_types_distribution(filtered_data())
    })
    output$mutation_types_distribution <- renderPlot({
      mutation_types_distribution_plot()
    })
    mutation_subtypes_distribution_plot <- reactive({
      mutation_subtypes_distribution(filtered_data())
    })
    output$mutation_subtypes_distribution <- renderPlot({
      mutation_subtypes_distribution_plot()
    })
    
    # SNP Analysis
    snp_class_comparison_plot <- reactive({
      snp_class_comparison(filtered_data())
    })
    output$snp_class_comparison <- renderPlot({
      snp_class_comparison_plot()
    })
    transversion_transitions_plot <- reactive({
      transversion_transitions(filtered_data())
    })
    output$transversion_transitions <- renderPlot({
      transversion_transitions_plot()
    })
    
    # INDEL Analysis
    insertion_deletions_plot <- reactive({
      insertion_deletions(filtered_data())
    })
    output$insertion_deletions <- renderPlot({
      insertion_deletions_plot()
    })
    indel_len_boxplot_plot <- reactive({
      indel_len_boxplot(filtered_data())
    })
    output$indel_len_boxplot <- renderPlot({
      indel_len_boxplot_plot()
    })
    
    # Quality Analysis
    quality_boxplot_plot <- reactive({
      quality_boxplot(filtered_data())
    })
    output$quality_boxplot <- renderPlot({
      quality_boxplot_plot()
    })
    
    # Allele Frequency Analysis
    frequency_ridges_plot <- reactive({
      frequency_ridges(filtered_data())
    })
    output$frequency_ridges <- renderPlot({
      frequency_ridges_plot()
    })
    
    # Read Depth Analysis
    read_depth_plot <- reactive({
      read_depth(filtered_data())
    })
    output$read_depth <- renderPlot({
      read_depth_plot()
    })
    
    # Venn diagram
    venn_diagram_plot <- eventReactive(input$generate_venn, {
      file1 <- input$file1_select
      file2 <- input$file2_select
      file3 <- input$file3_select
      req(file1 != "-", file2 != "-")
      if(file3 == "-"){
        venn_diagram(file1, file2, filtered_data())
      }else{
        venn_diagram_3files(file1, file2, file3, filtered_data())
      }
    })
    output$venn_diagram <- renderPlot({
      withProgress(message = "Generating Venn diagram...", value = 0, {
        venn_diagram_plot()
      })
    })
    
    # PDF REPORT
    output$download_vcf_comparison_pdf <- downloadHandler(
      filename ="vcfComparisonReport.pdf",
      content = function(file) {
        withProgress(message = "Generating VCF comparison report...", value = 0, {
          temp_report <- file.path("VCFComparisonModule", "vcfComparisonReport.Rmd")
          file.copy("vcfComparisonReport.Rmd", temp_report, overwrite = TRUE)
          incProgress(0.2, detail = "Loading comparison data...")
          venn_ready <- input$file1_select != "-" && input$file2_select != "-"
          params <- list(
            folder_name = vcf_comparison_folder_path(),
            selected_family = input$family_select,
            num_of_mutation = num_of_mutation_plot(),
            mutation_heatmap = mutation_heatmap_plot(),
            mutation_types_distribution = mutation_types_distribution_plot(),
            mutation_subtypes_distribution = mutation_subtypes_distribution_plot(),
            snp_class_comparison = snp_class_comparison_plot(),
            transversion_transitions = transversion_transitions_plot(),
            insertion_deletions = insertion_deletions_plot(),
            indel_len_boxplot = indel_len_boxplot_plot(),
            quality_boxplot = quality_boxplot_plot(),
            frequency_ridges = frequency_ridges_plot(),
            read_depth = read_depth_plot(),
            file1 = input$file1_select,
            file2 = input$file2_select,
            file3 = input$file3_select,
            venn_diagram = if (venn_ready) venn_diagram_plot() else NA,
            venn_ready = venn_ready
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
            title = "Venn diagram shows the number of shared and unique mutations between two VCF files based on matching chromosome, position, reference, and alternate alleles.",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        tags$h5("Choose Files"),
        fluidRow(
          column(3, selectInput(inputId = ns("file1_select"), 
                                label = "Select File 1 - VCF1:", 
                                choices = c(""), 
                                selected = "")),
          column(3, selectInput(inputId = ns("file2_select"), 
                                label = "Select File 2 - VCF2:", 
                                choices = c(""), 
                                selected = "")),
          column(3, selectInput(inputId = ns("file3_select"), 
                                label = "Select File 3 - VCF3 (Optional):", 
                                choices = c(""), 
                                selected = "")),
          column(3, actionButton(ns("generate_venn"), "Generate Venn Diagram"))
        ),
        fluidRow(column(12, plotOutput(ns("venn_diagram")))),
        downloadButton(ns("download_vcf_comparison_pdf"), "Download PDF Report"),
      )
    })
    
  })
}

