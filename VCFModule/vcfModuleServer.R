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
    filtered_data <- reactive({
      data <- processed_data()
      req(data)
      if (input$chrom_select != "All") {data %<>% filter(CHROM == input$chrom_select)}
      data
    })
    
    # Mutation counts
    mutation_donut_plot <- reactive({
      mutation_donut(filtered_data(), input$analysis_level == "Subtypes")
    })
    output$mutation_donut <- renderPlot({
      mutation_donut_plot()
    })
    mutation_distribution_plot <- reactive({
      mutation_distribution(filtered_data(), input$analysis_level == "Subtypes")
    })
    output$mutation_distribution <- renderPlot({
      mutation_distribution_plot()
    })
    
    # SNP Analysis
    output$snp_value <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      snp_value(data)
    })
    snp_types_donut_plot <- reactive({
      snp_types_donut(filtered_data())
    })
    output$snp_types_donut <- renderPlot({
      snp_types_donut_plot()
    })
    snp_class_stacked_plot <- reactive({
      snp_class_stacked(filtered_data())
    })
    output$snp_class_stacked <- renderPlot({
      snp_class_stacked_plot()
    })
    snp_class_boxplot_plot <- reactive({
      snp_class_boxplot(filtered_data())
    })
    output$snp_class_boxplot <- renderPlot({
      snp_class_boxplot_plot()
    })
    snp_class_barplot_plot <- reactive({
      snp_class_barplot(filtered_data())
    })
    output$snp_class_barplot <- renderPlot({
      snp_class_barplot_plot()
    })
    
    # INDEL Analysis
    output$indel_values <- renderText({
      data <- processed_data()
      if(input$chrom_select != "All"){data %<>% filter(CHROM == input$chrom_select)}
      indel_values(data)
    })
    indel_types_plot <- reactive({
      indel_types(filtered_data())
    })
    output$indel_types <- renderPlot({
      indel_types_plot()
    })
    indel_stacked_plot <- reactive({
      indel_stacked(filtered_data())
    })
    output$indel_stacked <- renderPlot({
      indel_stacked_plot()
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
    indel_length_plot <- reactive({
      indel_length(filtered_data())
    })
    output$indel_length <- renderPlotly({
      ggplotly(indel_length_plot())
    })
    indel_length_boxplot_plot <- reactive({
      indel_length_boxplot(filtered_data())
    })
    output$indel_length_boxplot <- renderPlot({
      indel_length_boxplot_plot()
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
    quality_bar_plot <- reactive({
      quality_bar(filtered_data())
    })
    output$quality_bar <- renderPlotly({
      ggplotly(quality_bar_plot())
    })
    quality_on_chroms_plot <- reactive({
      quality_on_chroms(filtered_data())
    })
    output$quality_on_chroms <- renderPlot({
      quality_on_chroms_plot()
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
    allele_freq_hexbin_plot <- reactive({
      allele_freq_hexbin(filtered_data())
    })
    output$allele_freq_hexbin <- renderPlotly({
      ggplotly(allele_freq_hexbin_plot())
    })
    allele_freq_on_chroms_plot <- reactive({
      allele_freq_on_chroms(filtered_data())
    })
    output$allele_freq_on_chroms <- renderPlot({
      allele_freq_on_chroms_plot()
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
    read_depth_density_plot <- reactive({
      read_depth_density(filtered_data())
    })
    output$read_depth_density <- renderPlotly({
      ggplotly(read_depth_density_plot())
    })
    read_depth_on_chroms_plot <- reactive({
      read_depth_on_chroms(filtered_data())
    })
    output$read_depth_on_chroms <- renderPlot({
      read_depth_on_chroms_plot()
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
            mutation_donut = mutation_donut_plot(),
            mutation_distribution = mutation_distribution_plot(),
            snp_value = snp_value(filtered_data()),
            snp_types_donut = snp_types_donut_plot(),
            snp_class_stacked = snp_class_stacked_plot(),
            snp_class_boxplot = snp_class_boxplot_plot(),
            snp_class_barplot = snp_class_barplot_plot(),
            indel_values = indel_values(filtered_data()),
            indel_length_avg = indel_length_avg(filtered_data()),
            indel_length_med = indel_length_med(filtered_data()),
            indel_types = indel_types_plot(),
            indel_stacked = indel_stacked_plot(),
            indel_length = indel_length_plot(),
            indel_length_boxplot = indel_length_boxplot_plot(),
            quality_avg = quality_avg(filtered_data()),
            quality_med = quality_med(filtered_data()),
            quality_bar = quality_bar_plot(),
            quality_on_chroms = quality_on_chroms_plot(),
            allele_freq_avg = allele_freq_avg(filtered_data()),
            allele_freq_med = allele_freq_med(filtered_data()),
            allele_freq_hexbin = allele_freq_hexbin_plot(),
            allele_freq_on_chroms = allele_freq_on_chroms_plot(),
            read_depth_avg = read_depth_avg(filtered_data()),
            read_depth_med = read_depth_med(filtered_data()),
            read_depth_density = read_depth_density_plot(),
            read_depth_on_chroms = read_depth_on_chroms_plot()
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
