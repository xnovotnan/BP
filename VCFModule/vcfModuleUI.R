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
      column(6, selectInput(ns("analysisLevel"), 
                            "Select Analysis Level:",
                            choices = c("Types", "Subtypes"), 
                            selected = "Types")),
      column(6, selectInput(ns("valueType"), 
                            "Select Value Type:",
                            choices = c("Absolute", "Percentage"), 
                            selected = "Percentage"))
    ),
    h4("Mutation counts"),
    fluidRow(
      column(6, plotOutput(ns("mutation_donut"))),
      column(6, plotOutput(ns("mutation_distribution")))
    ),
    h4("Single Nucleotide Polymorphism (SNP) Analysis"),
    fluidRow(
      column(12, verbatimTextOutput(ns("snp_values"))),
      column(4, plotOutput(ns("snp_types"))),
      column(8, plotOutput(ns("snp_class_stacked"))),
    ),
    fluidRow(
      column(7, plotOutput(ns("snp_class_combined"))),
      column(5, plotOutput(ns("snp_class")))
    ),
    h4("Insertion and Deletion (INDEL) Analysis"),
    fluidRow(
      column(12, verbatimTextOutput(ns("indel_values"))),
      column(4, plotOutput(ns("indel_types"))),
      column(8, plotOutput(ns("indel_stacked"))),
      column(6, verbatimTextOutput(ns("indel_length_avg"))),
      column(6, verbatimTextOutput(ns("indel_length_med"))),
      column(8, plotOutput(ns("indel_length"))),
      column(4, plotOutput(ns("indel_length_boxplot")))
    ),
    h4("Quality Analysis"),
    fluidRow(
      column(6, verbatimTextOutput(ns("qual_avg"))),
      column(6, verbatimTextOutput(ns("qual_med"))),
      column(6, plotOutput(ns("quality_bar"))),
      column(6, plotOutput(ns("quality_on_chroms")))
    ),
    h4("Allele Frequency Analysis"),
    fluidRow(
      column(6, verbatimTextOutput(ns("allele_freq_avg"))),
      column(6, verbatimTextOutput(ns("allele_freq_med"))),
      column(6, plotOutput(ns("allele_freq_hexbin"))),
      column(6, plotOutput(ns("allele_freq_on_chroms")))
    ),
    h4("Read Depth Analysis"),
    fluidRow(
      column(6, verbatimTextOutput(ns("read_depth_avg"))),
      column(6, verbatimTextOutput(ns("read_depth_med"))),
      column(6, plotOutput(ns("read_depth_density"))),
      column(6, plotOutput(ns("read_depth_on_chroms")))
    )
    
  )
}
