library(shiny)
library(bslib)
source(file.path("QualimapModule", "qualimapModuleUI.R"))
source(file.path("VCFModule", "vcfModuleUI.R"))
source(file.path("QualimapComparisonModule", "comparisonModuleUI.R"))
source(file.path("VCFComparisonModule", "vcfComparisonUI.R"))

ui <- page(
  title = "Dashboard App",
  div(class = "p-5",  
    navset_card_underline(
      nav_panel("VCF File Analysis", vcfModuleUI("vcf")),
      nav_panel("VCF Comparison", vcfComparisonUI("vcfComparison")),
      nav_panel("Qualimap Analysis", qualimapModuleUI("qualimap")),
      nav_panel("Qualimap Comparison", comparisonModuleUI("comparison"))
    )
  )
)
