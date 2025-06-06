library(shiny)
library(bslib)
source(file.path("VCFModule", "vcfModuleUI.R"))
source(file.path("VCFComparisonModule", "vcfComparisonUI.R"))
source(file.path("QualimapModule", "qualimapModuleUI.R"))
source(file.path("QualimapComparisonModule", "qualimapComparisonUI.R"))
source(file.path("FastqcModule", "fastqcModuleUI.R"))

# The ui function defines the main user interface in the application. 
# It displays the modules user interfaces in individual nav_panels, connecting 
# them with their server logic via the IDs.   

ui <- page(
  title = "Dashboard App",
  div(class = "p-5",  
    navset_card_underline(
      nav_panel("VCF File Analysis", vcfModuleUI("vcf")),
      nav_panel("VCF Comparison", vcfComparisonUI("vcfComparison")),
      nav_panel("Qualimap Analysis", qualimapModuleUI("qualimap")),
      nav_panel("Qualimap Comparison", qualimapComparisonUI("qualimapComparison")),
      nav_panel("FastQC Analysis", fastqcModuleUI("fastqc"))
    )
  )
)
