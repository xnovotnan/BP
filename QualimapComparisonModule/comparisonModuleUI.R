library(shiny)
library(shinyWidgets)
library(shinyFiles)

# UI komponent pre QUALIMAP porovnanie
comparisonModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Qualimap File Comparison"),
    shinyDirButton(ns("comparisonFolder"), "Select folder with qualimap files", "Select folder"),
    verbatimTextOutput(ns("selectedComparisonPath")),
    uiOutput(ns("qualimapModuleCompare"))
  )
}