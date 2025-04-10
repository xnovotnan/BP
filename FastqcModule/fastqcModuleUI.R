library(shiny)
library(shinyWidgets)
library(shinyFiles)

# UI komponent pre FastQC analýzu
fastqcModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3(
      "FastQC Statistics Summary",
      tags$span(
        icon("info-circle"),
        title = "Provides a detailed quality check of raw sequencing data. It helps identify potential issues such as low base quality, read duplication, contamination, or sequencing errors",
        style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
      )
    ),
    shinyDirButton(ns("fastqc_folder"), "Select fastqc folder", "Select folder"),
    verbatimTextOutput(ns("selected_fastqc")),
    uiOutput(ns("fastqc_module_combined"))
  )
}
