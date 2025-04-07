library(shiny)
library(shinyWidgets)
library(shinyFiles)

# UI komponent pre QUALIMAP porovnanie
comparisonModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3(
      "Qualimap File Comparison",
      tags$span(
        icon("info-circle"),
        title = "Compare alignment quality metrics (from Qualimap reports) across multiple samples",
        style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
      )
    ),
    shinyDirButton(ns("comparison_folder"), "Select folder with qualimap files", "Select folder"),
    verbatimTextOutput(ns("selected_comparison_path")),
    uiOutput(ns("qualimap_module_compare"))
  )
}