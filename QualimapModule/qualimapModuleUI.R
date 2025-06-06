library(shiny)
library(shinyWidgets)
library(shinyFiles)

# The qualimapModuleUI function defines the user interface for the QUALIMAP 
# analysis module, allowing users to select a QUALIMAP folder.  
qualimapModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    tags$h3(
      "Qualimap Analysis",
      tags$span(
        icon("info-circle"),
        title = "Analyze sequencing quality metrics from one sample, such as read mapping, coverage, GC content, and insert size",
        style = "cursor: help; color: black; font-size: 12px; vertical-align: middle;"
      )
    ),
    shinyDirButton(ns("qualimap_folder"), "Select qualimap folder", "Select folder"),
    verbatimTextOutput(ns("selected_folder_path")),
    uiOutput(ns("qualimap_module_combined"))
  )
}
