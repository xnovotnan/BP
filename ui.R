library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)
source("css.R")
source(file.path("QualimapModule", "qualimapModuleUI.R"))
source(file.path("VCFModule", "vcfModuleUI.R"))


ui <- page_sidebar(
  # ---------- SIDE BAR ----------
  title = "Uploading Files",
  sidebar = sidebar(
    theme = "bootstrap5",
    
    fileInput(inputId = "vcfFile", label = "Select VCF File", multiple = FALSE, 
              accept = ".vcf"),
    selectInput(inputId = "chrom_select", label = "Select Chromosome:", 
                choices = c("All"), selected = "All"),
    numericInput(inputId = "obs", label = "Number of observations to view:", 
                 value = 10),
    
    tags$hr(),
    h6("Select Qualimap Folder"),
    shinyDirButton("qualimapFolder", "Browse...", "Choose Qualimap Folder"),
    verbatimTextOutput("selectedFolderPath"),
  ),
  
  navset_card_underline(
    # ---------- VCF FILE SUMMARY WINDOW ----------
    nav_panel(
      "VCF File Summary", 
      fluidRow(
        column(12, div(verbatimTextOutput("file_summary"), 
                       style = "text-align: center; font-size: 16px; 
                               padding: 10px; border-radius: 5px;")),
        column(12, plotOutput("basic_visualizations", height = "800px")),
        column(12, div(dataTableOutput("table_view"), 
                       style="text-align: center; font-size: 16px;"))
      )
    ),
    
    # ---------- MUTATION ANALYSIS WINDOW ----------
    nav_panel(
      "Mutation Analysis",
      vcfModuleUI("vcf")
    ),
    
    # ---------- QUALIMAP ANALYSIS ----------
    nav_panel(
      "Qualimap Analysis",
      qualimapModuleUI("qualimap")
    )
  )
)