library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)

# UI komponent pre VCF analýzu
vcfModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Variant Summary per file"),
    fileInput(inputId = "vcfFile", 
              label = "Select VCF File", 
              multiple = FALSE, 
              accept = ".vcf"),
    uiOutput(ns("vcfModuleAnalysis"))
  )
}
