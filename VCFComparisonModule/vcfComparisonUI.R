library(shiny)
library(shinyWidgets)
library(shinyFiles)

# UI komponent pre VCF porovnanie
vcfComparisonUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("VCF File Comparison"),
    shinyDirButton(ns("vcfComparisonFolder"), "Select folder with VCF files", "Select folder"),
    verbatimTextOutput(ns("selectedComparisonPath")),
    uiOutput(ns("vcfModuleCompare"))
  )
}