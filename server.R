library(shiny)
library(bslib)
library(DT)
library(shinyFiles)
library(tinytex)
source("data_processing.R")
source("mutation_analysis.R")
source(file.path("QualimapModule", "qualimapAnalysis.R"))
source(file.path("QualimapModule", "qualimapModuleServer.R"))
options(shiny.maxRequestSize = 2000 * 1024^2)

# -----------------------------------
# SERVER
# -----------------------------------
# serverova cast shiny aplikacie, spracovanie vcf suboru, qualimap folderu 

server <- function(input, output, session) {
  # ---------- SIDE PANEL ----------
  
  # SELECT VCF FILE
  processedData <- reactive({
    req(input$vcf_file)
    withProgress(message = "Preprocessing file...", value = 0, {
      prepare_data(input$vcf_file$datapath)
    })
  })
  
  # SELECT CHROMOSOME
  observe({
    data <- processedData()
    req(data)
    chroms <- unique(data$CHROM)
    chrom_choices <- c("All", as.character(chroms))
    updateSelectInput(session, "chrom_select", choices = chrom_choices, selected = "All")
  })
  
  # SELECT QUALIMAP FOLDER
  volumes <- c(Home = path.expand("~"), Desktop = "~/Desktop", Documents = "~/Documents")
  shinyDirChoose(input, "qualimapFolder", roots = volumes, session = session)
  qualimapFolderPath <- reactive({
    req(typeof(input$qualimapFolder) != "integer")
    parseDirPath(volumes, input$qualimapFolder)
  })
  output$selectedFolderPath <- renderText({
    folderPath <- qualimapFolderPath()
    req(folderPath)
    paste(folderPath)
  })
  
  
  
  # ---------- FILE SUMMARY WINDOW ----------
  # summary statistika atributov
  output$file_summary <- renderPrint({
    data <- processedData()
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
      data <- processedData()
      req(data)
      if (input$chrom_select != "All") {
        data %<>% filter(CHROM == input$chrom_select)
      }
      plot_summary(data)
    }
  })
  
  # tabulka zo spracovaneho vstupneho suboru
  output$table_view <- renderDataTable({
    data <- processedData()
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
      data <- processedData()
      req(data)
      if (input$chrom_select != "All") {
        data %<>% filter(CHROM == input$chrom_select)
      }
      mut_summary(data, input$analysis_level == "Subtypes")
    }
  })
  
  # Heatmapa mutacii
  output$mut_heatmap <- renderPlot({
    data <- processedData()
    req(data)
    if(input$chrom_select != "All"){
      data %<>% filter(CHROM == input$chrom_select)
    }
    
    mutation_heatmap(data)
  })
  
  # Mutacie na chromozomoch
  output$mut_dist <- renderPlot({
    data <- processedData()
    req(data)
    if(input$chrom_select != "All"){
      data %<>% filter(CHROM == input$chrom_select)
    }
    
    mut_dist(data, input$analysis_level == "Subtypes")
  })
  
  
  # ---------- QUALIMAP ANALYSIS ----------
  qualimapModuleServer("qualimap", qualimapFolderPath)
  
}

