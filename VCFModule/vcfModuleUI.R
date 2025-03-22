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
    h5("Filters"),
    fluidRow(
      column(6, selectInput(ns("analysisLevel"), "Select Analysis Level:",
                            choices = c("Types", "Subtypes"), selected = "Types")),
      column(6, selectInput(ns("valueType"), "Select Value Type:",
                            choices = c("Absolute", "Percentage"), selected = "Absolute"))
    ),
    h4("Mutation counts"),
    fluidRow(
      column(6, plotOutput(ns("mut_summary"))),
      column(6, plotOutput(ns("mut_dist")))
    )
  )
}
