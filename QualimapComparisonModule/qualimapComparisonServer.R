library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("QualimapComparisonModule", "qualimapComparison.R"))

# The comparisonModuleServer function defines the server logic for the QUALIMAP comparison module.
# It processes the selected folder with QUALIMAP outputs, generates analysis components, updates the UI with 
# the results. It also generates PDF report from the analysis components. 
# The server function handles file reading, data processing, and renders the dynamic
# components based on the uploaded file.

comparisonModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(TUTUTUUTUTUUT = "~/Documents/BP/data/whitespring.qualimap", Home = path.expand("~"))
    shinyDirChoose(input, "comparison_folder", roots = volumes, session = session)

    comparison_folder_path <- reactive({
      req(typeof(input$comparison_folder) != "integer")
      parseDirPath(volumes, input$comparison_folder)
    })
    processed_data <- reactive({
      folder <- comparison_folder_path()
      req(folder)
      withProgress(message = "Preprocessing files...", value = 0, {
        process_qualimap_folders(folder)
      })
    })
  
    output$selected_comparison_path <- renderText({
      folder_path <- comparison_folder_path()
      req(folder_path)
      paste("Selected folder: ", folder_path)
    })
    
    # Reference
    output$bases_comparison <- renderPlot({
      bases_contigs_comparison(processed_data(), "number_of_bases", "Number of Bases (bp)")
    })
    output$contig_comparison <- renderPlot({
      bases_contigs_comparison(processed_data(), "number_of_contigs", "Number of Contigs")
    })
    
    # Read Statistics
    output$reads_comparison <- renderPlot({
      reads_comparison(processed_data())
    })
    output$mapped_paired_reads_comparison <- renderPlot({
      mapped_paired_reads(processed_data())
    })
    output$mapped_paired_reads_singletons_comparison <- renderPlot({
      mapped_paired_reads_singletons(processed_data())
    })
    output$mapped_bases_comparison <- renderPlot({
      mapped_bases_comparison(processed_data())
    })
    output$duplicated_reads_comparison <- renderPlot({
      duplicated_reads(processed_data())
    })
    output$mapping_quality_comparison <- renderPlot({
      mapping_quality_comparison(processed_data())
    })
    
    # Insert Size
    output$insert_size_mean_comparison <- renderPlot({
      insert_size_mean_comparison(processed_data())
    })
    
    # Data Coverage
    output$coverage_mean_comparison <- renderPlot({
      coverage_mean_comparison(processed_data())
    })
    
    # ACTG Content
    output$gc_percentage_comparison <- renderPlot({
      gc_comparison(processed_data())
    })
    output$actg_content_comparison <- renderPlot({
      stacked_actg(processed_data())
    })
    
    # PDF REPORT
    output$download_comparison_pdf <- downloadHandler(
      filename = "QUALIMAPcomparisonReport.pdf",
      content = function(file) {
        withProgress(message = "Generating comparison report...", value = 0, {
          temp_report <- file.path("QualimapComparisonModule", "qualimapComparisonReport.Rmd") 
          file.copy("qualimapReport.Rmd", temp_report, overwrite = TRUE)
          incProgress(0.2, detail = "Loading comparison data...")
          params <- list(
            folder_name = comparison_folder_path(),
            bases_comparison = bases_contigs_comparison(processed_data(), "number_of_bases", "Number of Bases (bp)"),
            contig_comparison = bases_contigs_comparison(processed_data(), "number_of_contigs", "Number of Contigs"),
            reads_comparison = reads_comparison(processed_data()),
            mapped_paired_reads_comparison = mapped_paired_reads(processed_data()),
            mapped_paired_reads_singletons_comparison = mapped_paired_reads_singletons(processed_data()),
            mapped_bases_comparison = mapped_bases_comparison(processed_data()),
            duplicated_reads_comparison = duplicated_reads(processed_data()),
            mapping_quality_comparison = mapping_quality_comparison(processed_data()),
            insert_size_mean_comparison = insert_size_mean_comparison(processed_data()),
            coverage_mean_comparison = coverage_mean_comparison(processed_data()),
            gc_percentage_comparison = gc_comparison(processed_data()),
            actg_content_comparison = stacked_actg(processed_data())
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
    
    
    # COMBINED QUALIMAP COMPARISON
    output$qualimap_module_compare <- renderUI({
      req(processed_data())
      ns <- NS(id)
      tagList(
        tags$h4(
          "Reference",
          tags$span(
            icon("info-circle"),
            title = "Shows details about the reference genome used for alignment",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, plotOutput(ns("bases_comparison"))),
          column(6, plotOutput(ns("contig_comparison")))
        ),
        tags$hr(),
        tags$h4(
          "Read Statistics",
          tags$span(
            icon("info-circle"),
            title = "Displays overview of read counts, mapping quality, and duplication metrics",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(12, plotOutput(ns("reads_comparison"))),
          tags$hr(),
          column(6, plotOutput(ns("mapped_paired_reads_comparison"))),
          column(6, plotOutput(ns("mapped_paired_reads_singletons_comparison"))),
          tags$hr(),
          column(12, plotOutput(ns("mapped_bases_comparison"))),
          tags$hr(),
          column(6, plotOutput(ns("duplicated_reads_comparison"))),
          column(6, plotOutput(ns("mapping_quality_comparison"))),
        ),
        tags$hr(),
        tags$h4(
          "Insert Size",
          tags$span(
            icon("info-circle"),
            title = "Displays distribution of insert sizes between paired-end reads",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(12, plotOutput(ns("insert_size_mean_comparison")))
        ),
        tags$hr(),
        tags$h4(
          "Data Coverage",
          tags$span(
            icon("info-circle"),
            title = "Shows genome coverage statistics across samples",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(12, plotOutput(ns("coverage_mean_comparison")))
        ),
        tags$hr(),
        tags$h4(
          "ACTG Content",
          tags$span(
            icon("info-circle"),
            title = "Displays GC percentage and base composition (A, C, T, G)",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, plotOutput(ns("gc_percentage_comparison"))),
          column(6, plotOutput(ns("actg_content_comparison"))),
        ),
        downloadButton(ns("download_comparison_pdf"), "Download PDF Report")
      )
    })
  })
}

