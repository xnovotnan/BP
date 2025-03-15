library(shiny)
library(bslib)
library(DT)
library(shinyFiles)
source("data_processing.R")
source("mutation_analysis.R")
source("qualimap_analysis.R")
options(shiny.maxRequestSize = 2000 * 1024^2)

# -----------------------------------
# SERVER
# -----------------------------------
# serverova cast shiny aplikacie, spracovanie vcf suboru, 

server <- function(input, output, session) {
  
  # ---------- SPRACOVANIE VCF SUBORU ----------
  processed_data <- reactive({
    req(input$file1)
    withProgress(message = "Preprocessing file...", value = 0, {
      prepare_data(input$file1$datapath)
    })
  })
  
  # ---------- SIDE PANEL ----------
  # dynamicky vyber chromozomu v bocnom paneli
  observe({
    data <- processed_data()
    req(data)
    chroms <- unique(data$CHROM)
    chrom_choices <- c("All", as.character(chroms))
    updateSelectInput(session, "chrom_select", choices = chrom_choices, selected = "All")
  })
  
  # ---------- FILE SUMMARY WINDOW ----------
  # zakladna summary statistika
  output$file_summary <- renderPrint({
    data <- processed_data()
    if (input$chrom_select != "All") {
      data %<>% filter(CHROM == input$chrom_select)
    }
    summary(data)
  })
  
  # grafy distribucie 
  output$basic_visualizations <- renderPlot({
    if (is.null(input$file1)) {
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
  
  # ukazka dat zo spracovaneho vstupneho suboru
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
    data <- processed_data()
    req(data)
    if (input$chrom_select != "All") {
      data %<>% filter(CHROM == input$chrom_select)
    }
    
    mut_summary(data, input$analysis_level == "Subtypes")
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
  # Selected folder path
  output$folder_path <- renderText({
    req(typeof(input$qualimap_folder) != "integer")
    folder_path <- parseDirPath(volumes, input$qualimap_folder)
    paste(folder_path)
  })
  
  # ELEMENTARY METRICS
  output$num_of_bases <- renderText({
    paste("Number of Bases:", qualimap_data()["number of bases"])
  })
  output$num_of_contigs <- renderText({
    paste("Number of Contigs:", qualimap_data()["number of contigs"])
  })
  output$mean_mapping_quality <- renderText({
    paste("Mean Mapping Quality:", qualimap_data()["mean mapping quality"])
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
  output$num_of_mapped_bases <- renderText({
    paste("Number of Mapped Bases:", qualimap_data()["number of mapped bases"])
  })
  output$num_of_sequenced_bases <- renderText({
    paste("Number of Sequenced Bases:", qualimap_data()["number of sequenced bases"])
  })
  output$num_of_duplicated_reads <- renderText({
    paste("Number of Duplicated Reads:", qualimap_data()["number of duplicated reads (flagged)"])
  })
  
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
  
  # DATA COVERAGE ANALYSIS
  output$mean_coverage <- renderText({
    paste("Mean Coverage:", qualimap_data()["mean coverageData"])
  })
  output$std_coverage <- renderText({
    paste("Standard Deviation of Coverage:", qualimap_data()["std coverageData"])
  })
  # Coverage graf
  output$qualimap_coverage <- renderPlot({
    req(typeof(input$qualimap_folder) != "integer")
    folder_path <- parseDirPath(volumes, input$qualimap_folder)
    process_qualimap_coverage(folder_path)
  })
  # Coverage per contig graf
  output$qualimap_coverage_pc <- renderPlot({
    req(typeof(input$qualimap_folder) != "integer")
    folder_path <- parseDirPath(volumes, input$qualimap_folder)
    process_qualimap_coverage_pc(folder_path)
  })
  
  # ACTG CONTENT
  output$gc_percentage <- renderText({
    paste("GC Percentage:", qualimap_data()["GC percentage"])
  })
  output$actg_content <- renderPlot({
    data <- qualimap_data()
    req(data)
    process_ACTG_content(data)
  })
  

}

