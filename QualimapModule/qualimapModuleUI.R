library(shiny)
library(shinyWidgets)
library(shinyFiles)

# UI komponent pre Qualimap analýzu
qualimapModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Qualimap Statistics Summary"),
    shinyDirButton(ns("qualimapFolder"), "Select qualimap folder", "Select folder"),
    verbatimTextOutput(ns("selectedFolderPath")),
    uiOutput(ns("qualimapModuleCombined"))
  )
}
