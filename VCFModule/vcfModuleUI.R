library(shiny)
library(shinyWidgets)
library(shinyFiles)
options(shiny.maxRequestSize = 2000 * 1024^2)

vcfModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Variant Summary per file"),
    fileInput(inputId = ns("vcfFile"), 
              label = "Select VCF File", 
              multiple = FALSE, 
              accept = ".vcf"),
    uiOutput(ns("vcfModuleAnalysis"))
  )
}
