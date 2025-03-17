library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)
source("data_processing.R")
source("qualimap_analysis.R")
options(shiny.maxRequestSize = 2000 * 1024^2)

# Serverová časť modulu pre Qualimap analýzu
qualimapModuleServer <- function(id, qualimapFolderPath) {
  moduleServer(id, function(input, output, session) {
    qualimapData <- reactive({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      process_qualimap(folderPath)
    })
    

    # REFERENCE
    output$num_of_bases <- renderText({
      paste("Number of Bases:", qualimapData()["number of bases"])
    })
    output$num_of_contigs <- renderText({
      paste("Number of Contigs:", qualimapData()["number of contigs"])
    })
    
    
    # READ STATISTICS
    output$num_of_reads <- renderText({
      paste("Number of Reads:", qualimapData()["number of reads"])
    })
    output$num_of_mapped_reads <- renderText({
      paste("Number of Mapped Reads:", qualimapData()["number of mapped reads"])
    })
    output$num_of_mapped_paired_reads <- renderText({
      paste("Number of Mapped Paired Reads:", qualimapData()["number of mapped paired reads (both in pair)"])
    })
    output$num_of_mapped_paired_reads_singletons <- renderText({
      paste("Number of Mapped Paired Reads (singletons):", qualimapData()["number of mapped paired reads (singletons)"])
    })
    output$num_of_mapped_bases <- renderText({
      paste("Number of Mapped Bases:", qualimapData()["number of mapped bases"])
    })
    output$num_of_sequenced_bases <- renderText({
      paste("Number of Sequenced Bases:", qualimapData()["number of sequenced bases"])
    })
    output$num_of_duplicated_reads <- renderText({
      paste("Number of Duplicated Reads:", qualimapData()["number of duplicated reads (flagged)"])
    })
    output$duplication_rate_histogram <- renderImage({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      file_path <- find_png(folderPath, "genome_uniq_read_starts_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # MAPPING QUALITY
    output$mean_mapping_quality <- renderText({
      paste("Mean Mapping Quality:", qualimapData()["mean mapping quality"])
    })
    output$mapping_quality_histogram <- renderImage({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      file_path <- find_png(folderPath, "genome_mapping_quality_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # INSERT SIZE
    output$mean_insert_size <- renderText({
      paste("Mean Insert Size:", qualimapData()["mean insert size"])
    })
    output$median_insert_size <- renderText({
      paste("Median Insert Size:", qualimapData()["median insert size"])
    })
    output$std_insert_size <- renderText({
      paste("Standard Deviation of Insert Size:", qualimapData()["std insert size"])
    })
    output$insert_size_across_reference <- renderImage({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      file_path <- find_png(folderPath, "genome_insert_size_across_reference.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    output$insert_size_histogram <- renderImage({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      file_path <- find_png(folderPath, "genome_insert_size_histogram.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # DATA COVERAGE ANALYSIS
    output$mean_coverage <- renderText({
      paste("Mean Coverage:", qualimapData()["mean coverageData"])
    })
    output$std_coverage <- renderText({
      paste("Standard Deviation of Coverage:", qualimapData()["std coverageData"])
    })
    output$qualimap_coverage <- renderPlot({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      process_qualimap_coverage(folderPath)
    })
    output$qualimap_coverage_pc <- renderPlot({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      process_qualimap_coverage_pc(folderPath)
    })
    
    
    # ACTG CONTENT
    output$gc_percentage <- renderText({
      paste("GC Percentage:", qualimapData()["GC percentage"])
    })
    output$actg_content_barplot <- renderPlot({
      data <- qualimapData()
      req(data)
      process_ACTG_content(data)
    })
    output$cg_content_distribution <- renderImage({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      file_path <- find_png(folderPath,"genome_gc_content_per_window.png")
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
          num_of_bases = qualimapData()["number of bases"],
          num_of_contigs = qualimapData()["number of contigs"],
          num_of_reads = qualimapData()["number of reads"],
          num_of_mapped_reads = qualimapData()["number of mapped reads"],
          num_of_mapped_paired_reads = qualimapData()["number of mapped paired reads (both in pair)"],
          num_of_mapped_paired_reads_singletons = qualimapData()["number of mapped paired reads (singletons)"],
          num_of_mapped_bases = qualimapData()["number of mapped bases"],
          num_of_sequenced_bases = qualimapData()["number of sequenced bases"],
          num_of_duplicated_reads = qualimapData()["number of duplicated reads (flagged)"],
          mean_insert_size = qualimapData()["mean insert size"],
          median_insert_size = qualimapData()["median insert size"],
          std_insert_size = qualimapData()["std insert size"],
          mean_coverage = qualimapData()["mean coverageData"],
          std_coverage = qualimapData()["std coverageData"],
          gc_percentage = qualimapData()["GC percentage"],
          duplication_rate_histogram = find_png(qualimapFolderPath(), "genome_uniq_read_starts_histogram.png"),
          mapping_quality_histogram = find_png(qualimapFolderPath(), "genome_mapping_quality_histogram.png"),
          insert_size_across_reference = find_png(qualimapFolderPath(), "genome_insert_size_across_reference.png"),
          insert_size_histogram = find_png(qualimapFolderPath(), "genome_insert_size_histogram.png"),
          qualimap_coverage = process_qualimap_coverage(qualimapFolderPath()),
          qualimap_coverage_pc = process_qualimap_coverage_pc(qualimapFolderPath()),  
          actg_content_barplot = process_ACTG_content(qualimapData()), 
          cg_content_distribution = find_png(qualimapFolderPath(),"genome_gc_content_per_window.png")
        )
        
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
  })
}
