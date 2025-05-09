library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("QualimapModule", "qualimapAnalysis.R"))

# The qualimapModuleServer function defines the server logic for the QUALIMAP analysis module.
# It processes the selected QUALIMAP folder, generates analysis components, updates the UI with 
# the results. It also generates PDF report from the analysis components. 
# The server function handles file reading, data processing, and renders the dynamic
# components based on the uploaded file.

qualimapModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(Home = path.expand("~"), Desktop = "~/Desktop", Documents = "~/Documents")
    shinyDirChoose(input, "qualimap_folder", roots = volumes, session = session)
    
    qualimap_folder_path <- reactive({
      req(typeof(input$qualimap_folder) != "integer")
      parseDirPath(volumes, input$qualimap_folder)
    })
    qualimap_data <- reactive({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      withProgress(message = "Preprocessing files...", value = 0, {
        process_qualimap(folder_path)
      })
    })
    output$selected_folder_path <- renderText({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      paste("Selected folder: ",folder_path)
    })

    # Reference
    output$bam_file <- renderText({
      paste("BAM File:", qualimap_data()["bam file"])
    })
    output$num_of_bases <- renderText({
      paste("Number of Bases:", qualimap_data()["number of bases"])
    })
    output$num_of_contigs <- renderText({
      paste("Number of Contigs:", qualimap_data()["number of contigs"])
    })
    
    # Read Statistics
    output$num_of_reads <- renderText({
      paste("Number of Reads:", qualimap_data()["number of reads"])
    })
    output$num_of_mapped_reads <- renderText({
      paste("Number of Mapped Reads:", qualimap_data()["number of mapped reads"])
    })
    output$num_of_mapped_paired_reads <- renderText({
      paste("Number of Mapped Paired Reads:", qualimap_data()["number of mapped paired reads (both in pair)"])
    })
    output$num_of_mapped_paired_reads_singletons <- renderText({
      paste("Number of Mapped Paired Reads (singletons):", qualimap_data()["number of mapped paired reads (singletons)"])
    })
    output$num_of_mapped_bases <- renderText({
      paste("Number of Mapped Bases:", qualimap_data()["number of mapped bases"])
    })
    output$num_of_sequenced_bases <- renderText({
      paste("Number of Sequenced Bases:", qualimap_data()["number of sequenced bases"])
    })
    output$num_of_duplicated_reads <- renderText({
      paste("Number of Duplicated Reads:", qualimap_data()["number of duplicated reads (flagged)"])
    })
    output$duplication_rate_histogram <- renderImage({
      file_path <- find_qualimap_png(qualimap_folder_path(), "genome_uniq_read_starts_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    output$mean_mapping_quality <- renderText({
      paste("Mean Mapping Quality:", qualimap_data()["mean mapping quality"])
    })
    output$mapping_quality_histogram <- renderImage({
      file_path <- find_qualimap_png(qualimap_folder_path(), "genome_mapping_quality_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    # Insert Size
    output$mean_insert_size <- renderText({
      paste("Mean Insert Size:", qualimap_data()["mean insert size"])
    })
    output$median_insert_size <- renderText({
      paste("Median Insert Size:", qualimap_data()["median insert size"])
    })
    output$std_insert_size <- renderText({
      paste("Standard Deviation of Insert Size:", qualimap_data()["std insert size"])
    })
    output$insert_size_across_reference <- renderImage({
      file_path <- find_qualimap_png(qualimap_folder_path(), "genome_insert_size_across_reference.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    output$insert_size_histogram <- renderImage({
      file_path <- find_qualimap_png(qualimap_folder_path(), "genome_insert_size_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # Data Coverage
    output$mean_coverage <- renderText({
      paste("Mean Coverage:", qualimap_data()["mean coverageData"])
    })
    output$std_coverage <- renderText({
      paste("Standard Deviation of Coverage:", qualimap_data()["std coverageData"])
    })
    output$qualimap_coverage <- renderPlot({
      process_qualimap_coverage(qualimap_folder_path())
    })
    output$qualimap_coverage_pc <- renderPlot({
      process_qualimap_coverage_pc(qualimap_folder_path())
    })
    
    
    # ACTG Content
    output$gc_percentage <- renderText({
      paste("GC Percentage:", qualimap_data()["GC percentage"])
    })
    output$actg_content_barplot <- renderPlot({
      process_ACTG_content(qualimap_data())
    })
    output$cg_content_distribution <- renderImage({
      file_path <- find_qualimap_png(qualimap_folder_path(),"genome_gc_content_per_window.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    

    # PDF REPORT
    output$download_qualimap_pdf <- downloadHandler(
      filename = "qualimapReport.pdf",
      content = function(file) {
        withProgress(message = "Generating Qualimap report...", value = 0, {
          temp_report <- file.path("QualimapModule", "qualimapReport.Rmd") 
          file.copy("qualimapReport.Rmd", temp_report, overwrite = TRUE)
          incProgress(0.2, detail = "Loading data...")
          params <- list(
            file_name = qualimap_folder_path(),
            bam_file = qualimap_data()["bam file"],
            num_of_bases = qualimap_data()["number of bases"],
            num_of_contigs = qualimap_data()["number of contigs"],
            num_of_reads = qualimap_data()["number of reads"],
            num_of_mapped_reads = qualimap_data()["number of mapped reads"],
            num_of_mapped_paired_reads = qualimap_data()["number of mapped paired reads (both in pair)"],
            num_of_mapped_paired_reads_singletons = qualimap_data()["number of mapped paired reads (singletons)"],
            num_of_mapped_bases = qualimap_data()["number of mapped bases"],
            num_of_sequenced_bases = qualimap_data()["number of sequenced bases"],
            num_of_duplicated_reads = qualimap_data()["number of duplicated reads (flagged)"],
            mean_insert_size = qualimap_data()["mean insert size"],
            median_insert_size = qualimap_data()["median insert size"],
            std_insert_size = qualimap_data()["std insert size"],
            mean_coverage = qualimap_data()["mean coverageData"],
            std_coverage = qualimap_data()["std coverageData"],
            gc_percentage = qualimap_data()["GC percentage"],
            duplication_rate_histogram = find_qualimap_png(qualimap_folder_path(), "genome_uniq_read_starts_histogram.png"),
            mapping_quality_histogram = find_qualimap_png(qualimap_folder_path(), "genome_mapping_quality_histogram.png"),
            insert_size_across_reference = find_qualimap_png(qualimap_folder_path(), "genome_insert_size_across_reference.png"),
            insert_size_histogram = find_qualimap_png(qualimap_folder_path(), "genome_insert_size_histogram.png"),
            qualimap_coverage = process_qualimap_coverage(qualimap_folder_path()),
            qualimap_coverage_pc = process_qualimap_coverage_pc(qualimap_folder_path()),  
            actg_content_barplot = process_ACTG_content(qualimap_data()), 
            cg_content_distribution = find_qualimap_png(qualimap_folder_path(), "genome_gc_content_per_window.png")
          )
          incProgress(0.5, detail = "Preparing report...")
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
    
    
    # COMBINED QUALIMAP ANALYSIS
    output$qualimap_module_combined <- renderUI({
      req(qualimap_folder_path())
      req(qualimap_data())
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
          column(12, verbatimTextOutput(ns("bam_file"))),
          column(6, verbatimTextOutput(ns("num_of_bases"))),
          column(6, verbatimTextOutput(ns("num_of_contigs")))
        ),
        tags$h4(
          "Read Statistics",
          tags$span(
            icon("info-circle"),
            title = "Displays overview of read counts, mapping quality, and duplication metrics",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, verbatimTextOutput(ns("num_of_reads"))),
          column(6, verbatimTextOutput(ns("num_of_mapped_reads"))),
          column(6, verbatimTextOutput(ns("num_of_mapped_paired_reads"))),
          column(6, verbatimTextOutput(ns("num_of_mapped_paired_reads_singletons"))),
          column(6, verbatimTextOutput(ns("num_of_mapped_bases"))),
          column(6, verbatimTextOutput(ns("num_of_sequenced_bases"))),
          column(6, verbatimTextOutput(ns("num_of_duplicated_reads"))),
          column(6, verbatimTextOutput(ns("mean_mapping_quality"))),
          column(6, imageOutput(ns("duplication_rate_histogram"), height = "500px")),
          column(6, imageOutput(ns("mapping_quality_histogram"), height = "500px"))
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
          column(6, verbatimTextOutput(ns("mean_insert_size"))),
          column(6, verbatimTextOutput(ns("median_insert_size"))),
          column(12, verbatimTextOutput(ns("std_insert_size"))),
          column(6, imageOutput(ns("insert_size_across_reference"), height = "500px")),
          column(6, imageOutput(ns("insert_size_histogram"), height = "500px"))
          
        ),
        tags$hr(),
        tags$h4(
          "Data Coverage",
          tags$span(
            icon("info-circle"),
            title = "Shows genome coverage statistics across sample",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, verbatimTextOutput(ns("mean_coverage"))),
          column(6, verbatimTextOutput(ns("std_coverage"))),
          column(6, plotOutput(ns("qualimap_coverage"))),
          column(6, plotOutput(ns("qualimap_coverage_pc")))
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
          column(12, verbatimTextOutput(ns("gc_percentage"))),
          column(6, plotOutput(ns("actg_content_barplot"), height = "500px")),
          column(6, imageOutput(ns("cg_content_distribution"), height = "500px"))
        ),
        tags$hr(),
        downloadButton(ns("download_qualimap_pdf"), "Download PDF Report")
      )
    })
  })
}
