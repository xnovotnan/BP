library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)
source("css.R")
source("qualimapModuleUI.R") 

ui <- page_sidebar(
  # ---------- SIDE BAR ----------
  title = "Uploading Files",
  sidebar = sidebar(
    theme = "bootstrap5",
    
    # Nacitanie suboru
    fileInput(inputId = "vcf_file", label = "Select VCF File", multiple = FALSE, 
              accept = ".vcf"),
    
    # Vyber chromozomu 
    selectInput(inputId = "chrom_select", label = "Select Chromosome:", 
                choices = c("All"), selected = "All"),
    
    # Pocet observacii
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
      fluidRow(
        column(12, h3("Variant Summary per file / per chromosome", 
                      style = "text-align: center; margin-top: 20px;")),
        selectInput(inputId = "analysis_level", label = "Select Analysis Level:",
                    choices = c("Types", "Subtypes"), selected = "Types")
      ),
      fluidRow(
        column(5, plotOutput("mut_summary")), 
        column(5, plotOutput("mut_dist"))
      ),
      fluidRow(
        column(5, plotOutput("mut_heatmap"))
      )
    ),
    
    
    # ---------- QUALIMAP ANALYSIS ----------
    nav_panel(
      "Qualimap Analysis",
      qualimapModuleUI("qualimap")
    )
  )
)