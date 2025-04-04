library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("ComparisonModule", "sampleComparison.R"))

# Serverová časť modulu pre Qualimap porovnanie
comparisonModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(TUTUTUUTUTUUT = "~/Documents/BP/data/whitespring.qualimap", Home = path.expand("~"))
    shinyDirChoose(input, "comparisonFolder", roots = volumes, session = session)

    comparisonFolderPath <- reactive({
      req(typeof(input$comparisonFolder) != "integer")
      parseDirPath(volumes, input$comparisonFolder)
    })
    processedData <- reactive({
      folder <- comparisonFolderPath()
      req(folder)
      process_qualimap_folders(folder)
    })
  
    output$selectedComparisonPath <- renderText({
      folderPath <- comparisonFolderPath()
      req(folderPath)
      paste("Selected folder: ",folderPath)
    })
    
    # Reference
    output$bases_comparison <- renderPlot({
      data <- processedData()
      req(data)
      bases_contigs_comparison(data, "number_of_bases", "Number of Bases (bp)")
    })
    output$contig_comparison <- renderPlot({
      data <- processedData()
      req(data)
      bases_contigs_comparison(data, "number_of_contigs", "Number of Contigs")
    })
    
    # Read Statistics
    output$reads_comparison <- renderPlot({
      data <- processedData()
      req(data)
      reads_comparison(data)
    })
    output$mapped_paired_reads_comparison <- renderPlot({
      data <- processedData()
      req(data)
      mapped_paired_reads(data)
    })
    output$mapped_paired_reads_singletons_comparison <- renderPlot({
      data <- processedData()
      req(data)
      mapped_paired_reads_singletons(data)
    })
    output$mapped_bases_comparison <- renderPlot({
      data <- processedData()
      req(data)
      mapped_bases_comparison(data)
    })
    output$duplicated_reads_comparison <- renderPlot({
      data <- processedData()
      req(data)
      duplicated_reads(data)
    })
    output$mapping_quality_comparison <- renderPlot({
      data <- processedData()
      req(data)
      mapping_quality_comparison(data)
    })
    
    # Insert Size
    output$insert_size_comparison <- renderPlot({
      data <- processedData()
      req(data)
      insert_size_comparison(data)
    })
    
    # Data Coverage
    output$coverage_comparison <- renderPlot({
      data <- processedData()
      req(data)
      coverage_comparison(data)
    })
    
    # ACTG Content
    output$gc_percentage_comparison <- renderPlot({
      data <- processedData()
      req(data)
      gc_comparison(data)
    })
    output$actg_content_comparison <- renderPlot({
      data <- processedData()
      req(data)
      stacked_actg(data)
    })
    
    # PDF REPORT
    output$download_comparison_pdf <- downloadHandler(
      filename = "comparisonReport.pdf",
      content = function(file) {
        tempReport <- file.path("ComparisonModule", "comparisonReport.Rmd") 
        file.copy("comparisonReport.Rmd", tempReport, overwrite = TRUE)
        params <- list(
          bases_comparison = bases_contigs_comparison(processedData(), "number_of_bases", "Number of Bases (bp)"),
          contig_comparison = bases_contigs_comparison(processedData(), "number_of_contigs", "Number of Contigs"),
          reads_comparison = reads_comparison(processedData()),
          mapped_paired_reads_comparison = mapped_paired_reads(processedData()),
          mapped_paired_reads_singletons_comparison = mapped_paired_reads_singletons(processedData()),
          mapped_bases_comparison = mapped_bases_comparison(processedData()),
          duplicated_reads_comparison = duplicated_reads(processedData()),
          mapping_quality_comparison = mapping_quality_comparison(processedData()),
          insert_size_comparison = insert_size_comparison(processedData()),
          coverage_comparison = coverage_comparison(processedData()),
          gc_percentage_comparison = gc_comparison(processedData()),
          actg_content_comparison = stacked_actg(processedData())
        )
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv()))
      }
    )
    
    # COMBINED QUALIMAP COMPARISON
    output$qualimapModuleCompare <- renderUI({
      req(processedData())
      ns <- NS(id)
      tagList(
        h4("Reference"),
        fluidRow(
          column(6, plotOutput(ns("bases_comparison"))),
          column(6, plotOutput(ns("contig_comparison")))
        ),
        tags$hr(),
        h4("Read Statistics"),
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
        h4("Insert Size"),
        fluidRow(
          column(12, plotOutput(ns("insert_size_comparison"))),
        ),
        tags$hr(),
        h4("Data Coverage"),
        fluidRow(
          column(12, plotOutput(ns("coverage_comparison")))
        ),
        tags$hr(),
        h4("ACTG Content"),
        fluidRow(
          column(6, plotOutput(ns("gc_percentage_comparison"))),
          column(6, plotOutput(ns("actg_content_comparison"))),
        ),
        downloadButton(ns("download_comparison_pdf"), "Download PDF Report")
      )
    })
  })
}

