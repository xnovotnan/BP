library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("QualimapModule", "qualimapAnalysis.R"))

# Serverová časť modulu pre Qualimap analýzu
qualimapModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(Home = path.expand("~"), Desktop = "~/Desktop", Documents = "~/Documents")
    shinyDirChoose(input, "qualimapFolder", roots = volumes, session = session)
    
    qualimapFolderPath <- reactive({
      req(typeof(input$qualimapFolder) != "integer")
      parseDirPath(volumes, input$qualimapFolder)
    })
    qualimapData <- reactive({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      process_qualimap(folderPath)
    })
    
    output$selectedFolderPath <- renderText({
      folderPath <- qualimapFolderPath()
      req(folderPath)
      paste("Selected folder: ",folderPath)
    })

    # Reference
    output$num_of_bases <- renderText({
      paste("Number of Bases:", qualimapData()["number of bases"])
    })
    output$num_of_contigs <- renderText({
      paste("Number of Contigs:", qualimapData()["number of contigs"])
    })
    
    
    # Read Statistics
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
    
    
    # Insert Size
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
    
    
    # Data Coverage
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
    
    
    # ACTG Content
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
      filename = "qualimapReport.pdf",
      content = function(file) {
        tempReport <- file.path("QualimapModule", "qualimapReport.Rmd") 
        file.copy("qualimapReport.Rmd", tempReport, overwrite = TRUE)
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
    
    
    # COMBINED QUALIMAP ANALYSIS
    output$qualimapModuleCombined <- renderUI({
      req(qualimapFolderPath())
      ns <- NS(id)
      tagList(
        h4("Reference"),
        fluidRow(
          column(6, verbatimTextOutput(ns("num_of_bases"))),
          column(6, verbatimTextOutput(ns("num_of_contigs")))
        ),
        h4("Read Statistics"),
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
        h4("Insert Size"),
        fluidRow(
          column(6, verbatimTextOutput(ns("mean_insert_size"))),
          column(6, verbatimTextOutput(ns("median_insert_size"))),
          column(12, verbatimTextOutput(ns("std_insert_size"))),
          column(6, imageOutput(ns("insert_size_across_reference"), height = "500px")),
          column(6, imageOutput(ns("insert_size_histogram"), height = "500px"))
          
        ),
        tags$hr(),
        h4("Data Coverage"),
        fluidRow(
          column(6, verbatimTextOutput(ns("mean_coverage"))),
          column(6, verbatimTextOutput(ns("std_coverage"))),
          column(6, plotOutput(ns("qualimap_coverage"))),
          column(6, plotOutput(ns("qualimap_coverage_pc")))
        ),
        tags$hr(),
        h4("ACTG Content"),
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
