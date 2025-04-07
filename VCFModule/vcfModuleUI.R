library(shiny)
library(shinyWidgets)
library(shinyFiles)
options(shiny.maxRequestSize = 2000 * 1024^2)

vcfModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3(
      "Variant Summary per file",
      tags$span(
        icon("info-circle"),
        title = "Explore detailed variant data from a single VCF file, including mutation types, quality, allele frequency, and read depth metrics",
        style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
      )
    ),
    fileInput(inputId = ns("vcf_file"), 
              label = "Select VCF File", 
              multiple = FALSE, 
              accept = ".gz"),
    uiOutput(ns("vcf_module_analysis"))
  )
}
