library(shiny)
library(shinyWidgets)
library(shinyFiles)
library(tinytex)
source(file.path("FastqcModule", "fastqcAnalysis.R"))

# Serverová časť modulu pre FastQC analýzu
fastqcModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    volumes <- c(TUTUTUUTUTUUT = "~/Documents/BP/data/sat.fastqc", Home = path.expand("~"), Desktop = "~/Desktop", Documents = "~/Documents")
    shinyDirChoose(input, "fastqc_folder", roots = volumes, session = session)
    
    fastqc_folder_path <- reactive({
      req(typeof(input$fastqc_folder) != "integer")
      parseDirPath(volumes, input$fastqc_folder)
    })
    fastqc_data_R1 <- reactive({
      folder_path <- fastqc_folder_path()
      req(folder_path)
      withProgress(message = "Preprocessing first read...", value = 0, {
        process_fastqc(folder_path, "R1")
      })
    })
    fastqc_data_R2 <- reactive({
      folder_path <- fastqc_folder_path()
      req(folder_path)
      withProgress(message = "Preprocessing second read...", value = 0, {
        process_fastqc(folder_path, "R2")
      })
    })
    
    output$selected_fastqc <- renderText({
      folder_path <- fastqc_folder_path()
      req(folder_path)
      paste0("Selected folder: ",folder_path)
    })
  
    # READ R1
    # Reference
    output$read_R1 <- renderText({
      paste0("Read: R1")
    })
    output$filename_R1 <- renderText({
      paste0("Filename: ", fastqc_data_R1()["Filename"])
    })
    output$file_type_R1 <- renderText({
      paste0("File type: ", fastqc_data_R1()["File type"])
    })
    output$encoding_R1 <- renderText({
      paste0("Encoding: ", fastqc_data_R1()["Encoding"])
    })
    output$total_sequences_R1 <- renderText({
      paste0("Total Sequences: ", fastqc_data_R1()["Total Sequences"])
    })
    output$total_bases_R1 <- renderText({
      paste0("Total Bases: ", fastqc_data_R1()["Total Bases"])
    })
    output$poor_quality_R1 <- renderText({
      paste0("Sequences flagged as poor quality: ", fastqc_data_R1()["Sequences flagged as poor quality"])
    })
    output$sequence_lenght_R1 <- renderText({
      paste0("Sequence length: ", fastqc_data_R1()["Sequence length"])
    })
    # Per base quality
    output$per_base_quality_R1 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R1", "per_base_quality.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Per sequence quality 
    output$per_sequence_quality_R1 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R1", "per_sequence_quality.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Per base sequence content
    output$per_base_sequenced_content_R1 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R1", "per_base_sequence_content.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Per base N content
    output$per_base_n_content_R1 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R1", "per_base_n_content.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Sequence Length Distribution
    output$sequence_length_distribution_R1 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R1", "sequence_length_distribution.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Sequence Duplication Levels
    output$duplication_levels_R1 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R1", "duplication_levels.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    
    
    # READ R2
    # Reference
    output$read_R2 <- renderText({
      paste0("Read: R2")
    })
    output$filename_R2 <- renderText({
      paste0("Filename: ", fastqc_data_R2()["Filename"])
    })
    output$file_type_R2 <- renderText({
      paste0("File type: ", fastqc_data_R2()["File type"])
    })
    output$encoding_R2 <- renderText({
      paste0("Encoding: ", fastqc_data_R2()["Encoding"])
    })
    output$total_sequences_R2 <- renderText({
      paste0("Total Sequences: ", fastqc_data_R2()["Total Sequences"])
    })
    output$total_bases_R2 <- renderText({
      paste0("Total Bases: ", fastqc_data_R2()["Total Bases"])
    })
    output$poor_quality_R2 <- renderText({
      paste0("Sequences flagged as poor quality: ", fastqc_data_R2()["Sequences flagged as poor quality"])
    })
    output$sequence_lenght_R2 <- renderText({
      paste0("Sequence length: ", fastqc_data_R2()["Sequence length"])
    })
    # Per base quality
    output$per_base_quality_R2 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R2", "per_base_quality.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Per sequence quality 
    output$per_sequence_quality_R2 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R2", "per_sequence_quality.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Per base sequence content
    output$per_base_sequenced_content_R2 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R2", "per_base_sequence_content.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Per base N content
    output$per_base_n_content_R2 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R2", "per_base_n_content.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Sequence Length Distribution
    output$sequence_length_distribution_R2 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R2", "sequence_length_distribution.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    # Sequence Duplication Levels
    output$duplication_levels_R2 <- renderImage({
      file_path <- find_fastqc_png(fastqc_folder_path(), "R2", "duplication_levels.png")
      list(src = file_path,
           contentType = "image/png",
           width = "100%")
    }, deleteFile = FALSE)
    

    # PDF REPORT
    output$download_fastqc_pdf <- downloadHandler(
      filename = "fastqcReport.pdf",
      content = function(file) {
        withProgress(message = "Generating FASTQC report...", value = 0, {
          temp_report <- file.path("FastqcModule", "fastqcReport.Rmd")
          file.copy("fastqcReport.Rmd", temp_report, overwrite = TRUE)
          incProgress(0.2, detail = "Loading FASTQC data...")
          params <- list(
            filename = paste0(fastqc_folder_path()),
            read_R1 = "Read: R1",
            filename_R1 = fastqc_data_R1()["Filename"],
            file_type_R1 = fastqc_data_R1()["File type"],
            encoding_R1 = fastqc_data_R1()["Encoding"],
            total_sequences_R1 = fastqc_data_R1()["Total Sequences"],
            total_bases_R1 = fastqc_data_R1()["Total Bases"],
            poor_quality_R1 = fastqc_data_R1()["Sequences flagged as poor quality"],
            sequence_length_R1 = fastqc_data_R1()["Sequence length"],
            per_base_quality_R1 = find_fastqc_png(fastqc_folder_path(), "R1", "per_base_quality.png"),
            per_sequence_quality_R1 = find_fastqc_png(fastqc_folder_path(), "R1", "per_sequence_quality.png"),
            per_base_sequence_content_R1 = find_fastqc_png(fastqc_folder_path(), "R1", "per_base_sequence_content.png"),
            per_base_n_content_R1 = find_fastqc_png(fastqc_folder_path(), "R1", "per_base_n_content.png"),
            sequence_length_distribution_R1 = find_fastqc_png(fastqc_folder_path(), "R1", "sequence_length_distribution.png"),
            duplication_levels_R1 = find_fastqc_png(fastqc_folder_path(), "R1", "duplication_levels.png"),
            
            read_R2 = "Read: R2",
            filename_R2 = fastqc_data_R2()["Filename"],
            file_type_R2 = fastqc_data_R2()["File type"],
            encoding_R2 = fastqc_data_R2()["Encoding"],
            total_sequences_R2 = fastqc_data_R2()["Total Sequences"],
            total_bases_R2 = fastqc_data_R2()["Total Bases"],
            poor_quality_R2 = fastqc_data_R2()["Sequences flagged as poor quality"],
            sequence_length_R2 = fastqc_data_R2()["Sequence length"],
            per_base_quality_R2 = find_fastqc_png(fastqc_folder_path(), "R2", "per_base_quality.png"),
            per_sequence_quality_R2 = find_fastqc_png(fastqc_folder_path(), "R2", "per_sequence_quality.png"),
            per_base_sequence_content_R2 = find_fastqc_png(fastqc_folder_path(), "R2", "per_base_sequence_content.png"),
            per_base_n_content_R2 = find_fastqc_png(fastqc_folder_path(), "R2", "per_base_n_content.png"),
            sequence_length_distribution_R2 = find_fastqc_png(fastqc_folder_path(), "R2", "sequence_length_distribution.png"),
            duplication_levels_R2 = find_fastqc_png(fastqc_folder_path(), "R2", "duplication_levels.png")
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
    

    # COMBINED fastqc ANALYSIS
    output$fastqc_module_combined <- renderUI({
      req(fastqc_folder_path())
      req(fastqc_data_R1())
      ns <- NS(id)
      tagList(
        tags$h4(
          "Reference",
          tags$span(
            icon("info-circle"),
            title = "Displays general file information and overall sequencing statistics",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, verbatimTextOutput(ns("read_R1"))),
          column(6, verbatimTextOutput(ns("read_R2"))),
          column(6, verbatimTextOutput(ns("filename_R1"))),
          column(6, verbatimTextOutput(ns("filename_R2"))),
          column(6, verbatimTextOutput(ns("file_type_R1"))),
          column(6, verbatimTextOutput(ns("file_type_R2"))),
          column(6, verbatimTextOutput(ns("encoding_R1"))),
          column(6, verbatimTextOutput(ns("encoding_R2"))),
          column(6, verbatimTextOutput(ns("total_sequences_R1"))),
          column(6, verbatimTextOutput(ns("total_sequences_R2"))),
          column(6, verbatimTextOutput(ns("total_bases_R1"))),
          column(6, verbatimTextOutput(ns("total_bases_R2"))),
          column(6, verbatimTextOutput(ns("poor_quality_R1"))),
          column(6, verbatimTextOutput(ns("poor_quality_R2"))),
          column(6, verbatimTextOutput(ns("sequence_lenght_R1"))),
          column(6, verbatimTextOutput(ns("sequence_lenght_R2")))
        ),
        tags$hr(),
        tags$h4(
          "Per Base Quality",
          tags$span(
            icon("info-circle"),
            title = "Shows how the quality score changes across each base position",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, imageOutput(ns("per_base_quality_R1"), height = "500px")),
          column(6, imageOutput(ns("per_base_quality_R2"), height = "500px"))
        ),
        tags$hr(),
        tags$h4(
          "Per Sequence Quality",
          tags$span(
            icon("info-circle"),
            title = "Displays the distribution of average quality scores per read",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, imageOutput(ns("per_sequence_quality_R1"), height = "500px")),
          column(6, imageOutput(ns("per_sequence_quality_R2"), height = "500px"))
          ),
        tags$hr(),
        tags$h4(
          "Per Base Sequence Content",
          tags$span(
            icon("info-circle"),
            title = "Shows how the proportion of each nucleotide (A, T, G, C) varies by position",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, imageOutput(ns("per_base_sequenced_content_R1"), height = "500px")),
          column(6, imageOutput(ns("per_base_sequenced_content_R2"), height = "500px"))
        ),
        tags$hr(),
        tags$h4(
          "Per Base N Content",
          tags$span(
            icon("info-circle"),
            title = "Displays the frequency of uncalled bases (“N”) across positions",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, imageOutput(ns("per_base_n_content_R1"),height = "500px")),
          column(6, imageOutput(ns("per_base_n_content_R2"),height = "500px"))
          ),
        tags$hr(),
        tags$h4(
          "Sequence Length Distribution",
          tags$span(
            icon("info-circle"),
            title = "Shows how read lengths are distributed across the dataset",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, imageOutput(ns("sequence_length_distribution_R1"),height = "500px")),
          column(6, imageOutput(ns("sequence_length_distribution_R2"),height = "500px"))
          ),
        tags$hr(),
        tags$h4(
          "Duplication Levels",
          tags$span(
            icon("info-circle"),
            title = "Displays how often sequences are duplicated",
            style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
          )
        ),
        fluidRow(
          column(6, imageOutput(ns("duplication_levels_R1"),height = "500px")),
          column(6, imageOutput(ns("duplication_levels_R2"),height = "500px"))
          ),
        downloadButton(ns("download_fastqc_pdf"), "Download PDF Report")
      )
    })
  })
}
