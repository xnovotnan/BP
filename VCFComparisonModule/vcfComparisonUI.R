library(shiny)
library(shinyWidgets)
library(shinyFiles)

# UI komponent pre VCF porovnanie
vcfComparisonUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3(
      "VCF File Comparison",
      tags$span(
        icon("info-circle"),
        title = "Compare the distribution and classification of somatic mutations, quality, allele frequency, and read depth analysis of multiple files",
        style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
      )
    ),
    shinyDirButton(ns("vcf_comparison_folder"), "Select folder with VCF files", "Select folder"),
    verbatimTextOutput(ns("selected_comparison_path")),
    uiOutput(ns("vcf_module_compare"))
  )
}