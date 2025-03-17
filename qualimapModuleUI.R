library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)

# UI komponent pre Qualimap analýzu
qualimapModuleUI <- function(id) {
  ns <- NS(id)
  uiOutput("qualimapModule")
  tagList(
    h3("Qualimap Statistics Summary"),
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
      column(6, imageOutput(ns("duplication_rate_histogram"))),
      column(6, imageOutput(ns("mapping_quality_histogram")))
    ),
    h4("Insert Size"),
    fluidRow(
      column(6, verbatimTextOutput(ns("mean_insert_size"))),
      column(6, verbatimTextOutput(ns("median_insert_size"))),
      column(12, verbatimTextOutput(ns("std_insert_size"))),
      column(6, imageOutput(ns("insert_size_across_reference"))),
      column(6, imageOutput(ns("insert_size_histogram")))
    ),
    h4("Data Coverage"),
    fluidRow(
      column(6, verbatimTextOutput(ns("mean_coverage"))),
      column(6, verbatimTextOutput(ns("std_coverage"))),
      column(6, plotOutput(ns("qualimap_coverage"))),
      column(6, plotOutput(ns("qualimap_coverage_pc")))
    ),
    h4("ACTG Content"),
    fluidRow(
      column(12, verbatimTextOutput(ns("gc_percentage"))),
      column(6, plotOutput(ns("actg_content_barplot"))),
      column(6, imageOutput(ns("cg_content_distribution")))
    ),
    downloadButton(ns("download_qualimap_pdf"), "Download PDF Report")
  )
}
