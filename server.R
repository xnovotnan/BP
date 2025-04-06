library(shiny)
source(file.path("VCFModule", "vcfModuleServer.R"))
source(file.path("QualimapModule", "qualimapModuleServer.R"))
source(file.path("QualimapComparisonModule", "comparisonModuleServer.R"))
source(file.path("VCFComparisonModule", "vcfComparisonServer.R"))

server <- function(input, output, session) {
  vcfModuleServer("vcf")
  qualimapModuleServer("qualimap")
  comparisonModuleServer("comparison")
  vcfComparisonServer("vcfComparison")
}
