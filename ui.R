library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)
source("css.R")

ui <- page_sidebar(
  # ---------- SIDE BAR ----------
  title = "Uploading VCF Files",
  sidebar = sidebar(
    theme = "bootstrap5",
    
    # Nacitanie suboru
    fileInput(inputId = "vcf_file", label = "Choose VCF File", multiple = FALSE, 
              accept = ".vcf"),
    
    # Vyber chromozomu 
    selectInput(inputId = "chrom_select", label = "Select Chromosome:", 
                choices = c("All"), selected = "All"),
    
    # Pocet observacii
    numericInput(inputId = "obs", label = "Number of observations to view:", 
                 value = 10),
    
    tags$hr(),
    
    h6("Select Qualimap Folder"),
    shinyDirButton("qualimap_folder", "Browse...", "Choose Qualimap Folder"),
    verbatimTextOutput("selected_folder_path"),
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
  )
)