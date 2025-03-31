library(shiny)
library(DT)
library(bslib)
library(shinyWidgets)
library(shinyFiles)
source("css.R")
source(file.path("QualimapModule", "qualimapModuleUI.R"))
source(file.path("VCFModule", "vcfModuleUI.R"))
source(file.path("ComparisonModule", "comparisonModuleUI.R"))


ui <- page(
  title = "Dashboard App",
  
  navset_card_underline(
    # ---------- MUTATION ANALYSIS WINDOW ----------
    nav_panel(
      "Mutation Analysis",
      vcfModuleUI("vcf")
    ),
    
    # ---------- QUALIMAP ANALYSIS ----------
    nav_panel(
      "Qualimap Analysis",
      qualimapModuleUI("qualimap")
    ),
    
    # ---------- SAMPLE COMPARISON ----------
    nav_panel(
      "Sample Comparison",
      comparisonModuleUI("comparison")
    )
  )
)