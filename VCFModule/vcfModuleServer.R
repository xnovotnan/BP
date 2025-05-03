library(shiny)
library(shinyWidgets)
library(tinytex)
source(file.path("VCFModule", "mutationAnalysis.R"))
library(plotly)

vcfModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    processed_data <- reactive({
      req(input$vcf_file)
      withProgress(message = "Preprocessing file...", value = 0, {
        prepare_data(input$vcf_file$datapath)
      })
    })
    observe({
      data <- processed_data()
      req(data)
      chroms <- unique(data$CHROM)
      chrom_choices <- c("All", as.character(chroms))
      updateSelectInput(session, "chrom_select", choices = chrom_choices, selected = "All")
    })
    
    # Mutation counts
    output$mutation_donut <- renderPlot({
      data <- processed_data()
      if (input$chrom_select != "All") {data %<>% filter(CHROM == input$chrom_select)}
      mutation_donut(data, input$analysis_level == "Subtypes")
    })
    output$mutation_distribution <- renderPlot({
      mutation_distribution(processed_data(), input$analysis_level == "Subtypes")
    })
    
    # SNP Analysis
    output$snp_value <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_value(data)
    })
    output$snp_types_donut <- renderPlot({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_types_donut(data)
    })
    output$snp_class_stacked <- renderPlot({
      snp_class_stacked(processed_data())
    })
    output$snp_class_boxplot <- renderPlot({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_class_boxplot(data)
    })
    output$snp_class_barplot <- renderPlot({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_class_barplot(data)
    })
    
    # INDEL Analysis
    output$indel_values <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_values(data)
    })
    output$indel_types <- renderPlot({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_types(data)
    })
    output$indel_stacked <- renderPlot({
      indel_stacked(processed_data())
    })
    output$indel_length_avg <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_length_avg(data)
    })
    output$indel_length_med <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_length_med(data)
    })
    output$indel_length <- renderPlotly({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      ggplotly(indel_length(data))
    })
    output$indel_length_boxplot <- renderPlot({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_length_boxplot(data)
    })
    
    # Quality
    output$qual_avg <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      quality_avg(data)
    })
    output$qual_med <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      quality_med(data)
    })
    output$quality_bar <- renderPlotly({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      ggplotly(quality_bar(data))
    })
    output$quality_on_chroms <- renderPlot({
      quality_on_chroms(processed_data())
    })
    
    # Allele frequency
    output$allele_freq_avg <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      allele_freq_avg(data)
    })
    output$allele_freq_med <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      allele_freq_med(data)
    })
    output$allele_freq_hexbin <- renderPlotly({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      ggplotly(allele_freq_hexbin(data))
    })
    output$allele_freq_on_chroms <- renderPlot({
      allele_freq_on_chroms(processed_data())
    })
      
    # Read Depth
    output$read_depth_avg <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      read_depth_avg(data)
    })
    output$read_depth_med <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      read_depth_med(data)
    })
    output$read_depth_density <- renderPlotly({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      ggplotly(read_depth_density(data))
    })
    output$read_depth_on_chroms <- renderPlot({
      read_depth_on_chroms(processed_data())
    })
    
    # PDF REPORT
    output$download_vcf_pdf <- downloadHandler(
      filename = function() {
        "vcfReport.pdf"
      },
      content = function(file) {
        withProgress(message = "Generating report...", value = 0, {
          temp_report <- file.path("VCFModule", "vcfReport.Rmd") 
          file.copy("vcfReport.Rmd", temp_report, overwrite = TRUE)
          
          params <- list(
            file_name = input$vcf_file[[1]],
            selected_chromosome = input$chrom_select,
            mutation_donut = mutation_donut(processed_data()),
            mutation_distribution = mutation_distribution(processed_data()),
            snp_value = snp_value(processed_data()),
            snp_types_donut = snp_types_donut(processed_data()),
            snp_class_stacked = snp_class_stacked(processed_data()),
            snp_class_boxplot = snp_class_boxplot(processed_data()),
            snp_class_barplot = snp_class_barplot(processed_data()),
            indel_values = indel_values(processed_data()),
            indel_length_avg = indel_length_avg(processed_data()),
            indel_length_med = indel_length_med(processed_data()),
            indel_types = indel_types(processed_data()),
            indel_stacked = indel_stacked(processed_data()),
            indel_length = indel_length(processed_data()),
            indel_length_boxplot = indel_length_boxplot(processed_data()),
            quality_avg = quality_avg(processed_data()),
            quality_med = quality_med(processed_data()),
            quality_bar = quality_bar(processed_data()),
            quality_on_chroms = quality_on_chroms(processed_data()),
            allele_freq_avg = allele_freq_avg(processed_data()),
            allele_freq_med = allele_freq_med(processed_data()),
            allele_freq_hexbin = allele_freq_hexbin(processed_data()),
            allele_freq_on_chroms = allele_freq_on_chroms(processed_data()),
            read_depth_avg = read_depth_avg(processed_data()),
            read_depth_med = read_depth_med(processed_data()),
            read_depth_density = read_depth_density(processed_data()),
            read_depth_on_chroms = read_depth_on_chroms(processed_data())
          )
          
          incProgress(0.4, detail = "Preparing data...")
          
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
    
    
    # COMBINED VCF ANALYSIS
    output$vcf_module_analysis <- renderUI({
      req(processed_data())
      ns <- NS(id)
      tagList(
        tags$h5(
          "Filters",
          tags$span(
            icon("info-circle"),
            title = "Focus the analysis on a specific chromosome to see its specifics",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, selectInput(ns("analysis_level"), 
                                "Select Analysis Level:",
                                choices = c("Types", "Subtypes"), 
                                selected = "Types")),
          column(6, selectInput(inputId = ns("chrom_select"), 
                                label = "Select Chromosome:", 
                                choices = c("All"), 
                                selected = "All"))
        ),
        tags$h4(
          "Mutation Count Overview",
          tags$span(
            icon("info-circle"),
            title = "Summarizes mutation types and their distribution to detect mutation hotspots",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, plotOutput(ns("mutation_donut"))),
          column(6, plotOutput(ns("mutation_distribution")))
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
          column(12, verbatimTextOutput(ns("snp_value"))),
          column(4, plotOutput(ns("snp_types_donut"))),
          column(8, plotOutput(ns("snp_class_stacked"))),
          tags$hr(),
          column(7, plotOutput(ns("snp_class_boxplot"))),
          column(5, plotOutput(ns("snp_class_barplot")))
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
          column(12, verbatimTextOutput(ns("indel_values"))),
          column(4, plotOutput(ns("indel_types"))),
          column(8, plotOutput(ns("indel_stacked"))),
          column(6, verbatimTextOutput(ns("indel_length_avg"))),
          column(6, verbatimTextOutput(ns("indel_length_med"))),
          column(8, plotlyOutput(ns("indel_length"))),
          column(4, plotOutput(ns("indel_length_boxplot")))
        ),
        tags$hr(),
        tags$h4(
          "Quality Analysis",
          tags$span(
            icon("info-circle"),
            title = "Displays quality scores, including average, median, and distribution by chromosome",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, verbatimTextOutput(ns("qual_avg"))),
          column(6, verbatimTextOutput(ns("qual_med"))),
          column(6, plotlyOutput(ns("quality_bar"))),
          column(6, plotOutput(ns("quality_on_chroms")))
        ),
        tags$hr(),
        tags$h4(
          "Allele Frequency Analysis",
          tags$span(
            icon("info-circle"),
            title = "Displays the distribution of allele frequency and how it varies across chromosomes",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, verbatimTextOutput(ns("allele_freq_avg"))),
          column(6, verbatimTextOutput(ns("allele_freq_med"))),
          column(6, plotlyOutput(ns("allele_freq_hexbin"))),
          column(6, plotOutput(ns("allele_freq_on_chroms")))
        ),
        tags$hr(),
        tags$h4(
          "Read Depth Analysis",
          tags$span(
            icon("info-circle"),
            title = "Displays the distribution of read depth and how it varies across chromosomes",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, verbatimTextOutput(ns("read_depth_avg"))),
          column(6, verbatimTextOutput(ns("read_depth_med"))),
          column(6, plotlyOutput(ns("read_depth_density"))),
          column(6, plotOutput(ns("read_depth_on_chroms")))
        ),
        downloadButton(ns("download_vcf_pdf"), "Download PDF Report")
      )
    })
  })
}
