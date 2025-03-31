library(shiny)
library(bslib)
library(DT)
library(shinyFiles)
library(tinytex)
source("data_processing.R")
source(file.path("VCFModule", "vcfModuleServer.R"))
source(file.path("QualimapModule", "qualimapModuleServer.R"))
source(file.path("ComparisonModule", "comparisonModuleServer.R"))
options(shiny.maxRequestSize = 2000 * 1024^2)

# -----------------------------------
# SERVER
# -----------------------------------
# serverova cast shiny aplikacie, spracovanie vcf suboru, qualimap folderu 

server <- function(input, output, session) {
  # ---------- SIDE PANEL ----------
  # SELECT VCF FILE
  processedData <- reactive({
    req(input$vcfFile)
    withProgress(message = "Preprocessing file...", value = 0, {
      prepare_data(input$vcfFile$datapath)
    })
  })
  
  # SELECT CHROMOSOME
  observe({
    data <- processedData()
    req(data)
    chroms <- unique(data$CHROM)
    chromChoices <- c("All", as.character(chroms))
    updateSelectInput(session, "chrom_select", choices = chromChoices, selected = "All")
  })
  
  chrom_select_val <- reactiveVal("All")
  observe({
    chrom_select_val(input$chrom_select)
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
    if (is.null(input$vcfFile)) {
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
  vcfModuleServer("vcf", processedData, chrom_select_val)
  
  # ---------- QUALIMAP ANALYSIS ----------
  qualimapModuleServer("qualimap")
  
  # ---------- SAMPLE COMPARISON ----------
  comparisonModuleServer("comparison")
  
}

