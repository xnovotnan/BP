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
    uiOutput(ns("qualimapModuleCombined"))
  )
}
