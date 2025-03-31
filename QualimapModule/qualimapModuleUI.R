library(shiny)
library(bslib)
library(shinyWidgets)
library(DT)
library(shinyFiles)
library(tinytex)

# UI komponent pre Qualimap analýzu
qualimapModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h6("Select Qualimap Folder"),
    shinyDirButton(ns("qualimapFolder"), "Browse...", "Choose Qualimap Folder"),
    verbatimTextOutput(ns("selectedFolderPath")),
    uiOutput(ns("qualimapModuleCombined"))
  )
}
