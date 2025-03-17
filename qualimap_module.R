library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)
source("data_processing.R")
source("mutation_analysis.R")
source("qualimap_analysis.R")
options(shiny.maxRequestSize = 2000 * 1024^2)

# UI komponent pre Qualimap analýzu
qualimapModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Qualimap Statistics Summary"),
    h4("Reference"),
    fluidRow(
      column(6, verbatimTextOutput("num_of_bases")),
      column(6, verbatimTextOutput("num_of_contigs"))
    ),
    h4("Read Statistics"),
    fluidRow(
      column(6, verbatimTextOutput("num_of_reads")),
      column(6, verbatimTextOutput("num_of_mapped_reads")),
      column(6, verbatimTextOutput("num_of_mapped_paired_reads")),
      column(6, verbatimTextOutput("num_of_mapped_paired_reads_singletons")),
      column(6, verbatimTextOutput("num_of_mapped_bases")),
      column(6, verbatimTextOutput("num_of_sequenced_bases")),
      column(6, verbatimTextOutput("num_of_duplicated_reads")),
      column(6, verbatimTextOutput("mean_mapping_quality")),
      column(6, imageOutput("duplication_rate_histogram")),
      column(6, imageOutput("mapping_quality_histogram"))
    ),
    h4("Insert Size"),
    fluidRow(
      column(6, verbatimTextOutput("mean_insert_size")),
      column(6, verbatimTextOutput("median_insert_size")),
      column(12, verbatimTextOutput("std_insert_size")),
      column(6, imageOutput("insert_size_across_reference")),
      column(6, imageOutput("insert_size_histogram"))
    ),
    h4("Data Coverage"),
    fluidRow(
      column(6, verbatimTextOutput("mean_coverage")),
      column(6, verbatimTextOutput("std_coverage")),
      column(6, plotOutput("qualimap_coverage")),
      column(6, plotOutput("qualimap_coverage_pc"))
    ),
    h4("ACTG Content"),
    fluidRow(
      column(12, verbatimTextOutput("gc_percentage")),
      column(6, plotOutput("actg_content_barplot")),
      column(6, imageOutput("cg_content_distribution"))
    ),
    downloadButton("download_qualimap_pdf", "Download PDF Report")
  )
}


# Serverová časť modulu pre Qualimap analýzu
qualimapModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Spracovanie Qualimap foldera/suborov
    volumes <- c(Home = path.expand("~"), Desktop = "~/Desktop", Documents = "~/Documents")
    shinyDirChoose(input, "qualimap_folder", roots = volumes, session = session)
    qualimap_data <- reactive({
      req(typeof(input$qualimap_folder) != "integer")
      folder_path <- parseDirPath(volumes, input$qualimap_folder)
      process_qualimap(folder_path)
    })
    qualimap_folder_path <- reactive({
      req(typeof(input$qualimap_folder) != "integer")
      folder_path <- parseDirPath(volumes, input$qualimap_folder)
    })
    output$selected_folder_path <- renderText({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      paste(folder_path)
    })
    
    
    # REFERENCE
    output$num_of_bases <- renderText({
      paste("Number of Bases:", qualimap_data()["number of bases"])
    })
    output$num_of_contigs <- renderText({
      paste("Number of Contigs:", qualimap_data()["number of contigs"])
    })
    
    
    # READ STATISTICS
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
      folder_path <- qualimap_folder_path()
      req(folder_path)
      file_path <- find_png(folder_path, "genome_uniq_read_starts_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # MAPPING QUALITY
    output$mean_mapping_quality <- renderText({
      paste("Mean Mapping Quality:", qualimap_data()["mean mapping quality"])
    })
    output$mapping_quality_histogram <- renderImage({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      file_path <- find_png(folder_path, "genome_mapping_quality_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # INSERT SIZE
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
      folder_path <- qualimap_folder_path()
      req(folder_path)
      file_path <- find_png(folder_path, "genome_insert_size_across_reference.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    output$insert_size_histogram <- renderImage({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      file_path <- find_png(folder_path, "genome_insert_size_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # DATA COVERAGE ANALYSIS
    output$mean_coverage <- renderText({
      paste("Mean Coverage:", qualimap_data()["mean coverageData"])
    })
    output$std_coverage <- renderText({
      paste("Standard Deviation of Coverage:", qualimap_data()["std coverageData"])
    })
    output$qualimap_coverage <- renderPlot({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      process_qualimap_coverage(folder_path)
    })
    output$qualimap_coverage_pc <- renderPlot({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      process_qualimap_coverage_pc(folder_path)
    })
    
    
    # ACTG CONTENT
    output$gc_percentage <- renderText({
      paste("GC Percentage:", qualimap_data()["GC percentage"])
    })
    output$actg_content_barplot <- renderPlot({
      data <- qualimap_data()
      req(data)
      process_ACTG_content(data)
    })
    output$cg_content_distribution <- renderImage({
      folder_path <- qualimap_folder_path()
      req(folder_path)
      file_path <- find_png(folder_path,"genome_gc_content_per_window.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    

    # PDF REPORT
    output$download_qualimap_pdf <- downloadHandler(
      filename = "report.pdf",
      content = function(file) {
        tempReport <- file.path(tempdir(), "qualimap_report.Rmd")
        file.copy("qualimap_report.Rmd", tempReport, overwrite = TRUE)
        params <- list(
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
          duplication_rate_histogram = find_png(qualimap_folder_path(), "genome_uniq_read_starts_histogram.png"),
          mapping_quality_histogram = find_png(qualimap_folder_path(), "genome_mapping_quality_histogram.png"),
          insert_size_across_reference = find_png(qualimap_folder_path(), "genome_insert_size_across_reference.png"),
          insert_size_histogram = find_png(qualimap_folder_path(), "genome_insert_size_histogram.png"),
          qualimap_coverage = process_qualimap_coverage(qualimap_folder_path()),
          qualimap_coverage_pc = process_qualimap_coverage_pc(qualimap_folder_path()),  
          actg_content_barplot = process_ACTG_content(qualimap_data()), 
          cg_content_distribution = find_png(qualimap_folder_path(),"genome_gc_content_per_window.png")
        )
        
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
  })
}
