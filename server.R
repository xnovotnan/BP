library(shiny)
library(bslib)
library(DT)
library(shinyFiles)
library(tinytex)
source("data_processing.R")
source("mutation_analysis.R")
source("qualimap_analysis.R")
options(shiny.maxRequestSize = 2000 * 1024^2)

# -----------------------------------
# SERVER
# -----------------------------------
# serverova cast shiny aplikacie, spracovanie vcf suboru, 

server <- function(input, output, session) {
  # ---------- SIDE PANEL ----------
  
  # SELECT VCF FILE
  processed_data <- reactive({
    req(input$vcf_file)
    withProgress(message = "Preprocessing file...", value = 0, {
      prepare_data(input$vcf_file$datapath)
    })
  })
  
  # SELECT CHROMOSOME
  observe({
    data <- processed_data()
    req(data)
    chroms <- unique(data$CHROM)
    chrom_choices <- c("All", as.character(chroms))
    updateSelectInput(session, "chrom_select", choices = chrom_choices, selected = "All")
  })
  
  
  # ---------- FILE SUMMARY WINDOW ----------
  # summary statistika atributov
  output$file_summary <- renderPrint({
    data <- processed_data()
    if (input$chrom_select != "All") {
      data %<>% filter(CHROM == input$chrom_select)
    }
    summary(data)
  })
  
  # grafy distribucie 
  output$basic_visualizations <- renderPlot({
    if (is.null(input$vcf_file)) {
      ggplot() + 
        theme_minimal() +
        ggtitle("Waiting for file upload...") + 
        theme(
          plot.title = element_text(size = 16, hjust = 0.5)
        )
    } else {
      data <- processed_data()
      req(data)
      if (input$chrom_select != "All") {
        data %<>% filter(CHROM == input$chrom_select)
      }
      plot_summary(data)
    }
  })
  
  # tabulka zo spracovaneho vstupneho suboru
  output$table_view <- renderDataTable({
    data <- processed_data()
    if (input$chrom_select != "All") {
      data %<>% filter(CHROM == input$chrom_select)
    }
    head(data, n = input$obs)
  }, options = list( 
    autoWidth = TRUE,
    columnDefs = list(list(targets = "_all", className = "dt-center small-font")),
    scrollX = FALSE,
    searching = FALSE,
    ordering = FALSE
  ),selection = "none")
  
  
  # ---------- MUTATION ANALYSIS WINDOW ----------
  # Zhrnutie mutacii
  output$mut_summary <- renderPlot({
    if (is.null(input$vcf_file)) {
      ggplot() + 
        theme_minimal() +
        ggtitle("Waiting for file upload...") + 
        theme(
          plot.title = element_text(size = 16, hjust = 0.5)
        )
    } else {
      data <- processed_data()
      req(data)
      if (input$chrom_select != "All") {
        data %<>% filter(CHROM == input$chrom_select)
      }
      mut_summary(data, input$analysis_level == "Subtypes")
    }
  })
  
  # Heatmapa mutacii
  output$mut_heatmap <- renderPlot({
    data <- processed_data()
    req(data)
    if(input$chrom_select != "All"){
      data %<>% filter(CHROM == input$chrom_select)
    }
    
    mutation_heatmap(data)
  })
  
  # Mutacie na chromozomoch
  output$mut_dist <- renderPlot({
    data <- processed_data()
    req(data)
    if(input$chrom_select != "All"){
      data %<>% filter(CHROM == input$chrom_select)
    }
    
    mut_dist(data, input$analysis_level == "Subtypes")
  })
  
  
  # ---------- QUALIMAP ANALYSIS ----------
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
  # Selected folder path
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
  # duplication rate histogram
  output$duplication_rate_histogram <- renderImage({
    folder_path <- qualimap_folder_path()
    req(folder_path)
    file_path <- show_duplication_rate_histogram(folder_path)
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
    file_path <- show_mapping_quality_histogram(folder_path)
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
  # insert size across reference
  output$insert_size_across_reference <- renderImage({
    folder_path <- qualimap_folder_path()
    req(folder_path)
    file_path <- show_insert_size_across_reference(folder_path)
    list(src = file_path,
         contentType = "image/png",
         width = "100%")
  }, deleteFile = FALSE)
  # insert size histogram
  output$insert_size_histogram <- renderImage({
    folder_path <- qualimap_folder_path()
    req(folder_path)
    file_path <- show_insert_size_histogram(folder_path)
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
  # Coverage graf
  output$qualimap_coverage <- renderPlot({
    folder_path <- qualimap_folder_path()
    req(folder_path)
    process_qualimap_coverage(folder_path)
  })
  # Coverage per contig graf
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
    file_path <- show_gc_content(folder_path)
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
        duplication_rate_histogram = show_duplication_rate_histogram(qualimap_folder_path()),
        mapping_quality_histogram = show_mapping_quality_histogram(qualimap_folder_path()),
        insert_size_across_reference = show_insert_size_across_reference(qualimap_folder_path()),
        insert_size_histogram = show_insert_size_histogram(qualimap_folder_path()),
        qualimap_coverage = process_qualimap_coverage(qualimap_folder_path()),
        qualimap_coverage_pc = process_qualimap_coverage_pc(qualimap_folder_path()),  
        actg_content_barplot = process_ACTG_content(qualimap_data()), 
        cg_content_distribution = show_gc_content(qualimap_folder_path())   
      )

      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
}

