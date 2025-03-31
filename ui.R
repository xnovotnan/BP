library(shiny)
library(bslib)
source(file.path("QualimapModule", "qualimapModuleUI.R"))
source(file.path("VCFModule", "vcfModuleUI.R"))
source(file.path("ComparisonModule", "comparisonModuleUI.R"))

ui <- page(
  title = "Dashboard App",
  div(class = "p-5",  
    navset_card_underline(
      nav_panel("Mutation Analysis", vcfModuleUI("vcf")),
      nav_panel("Qualimap Analysis", qualimapModuleUI("qualimap")),
      nav_panel("Sample Comparison", comparisonModuleUI("comparison"))
    )
  )
)
