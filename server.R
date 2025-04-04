library(shiny)
source(file.path("VCFModule", "vcfModuleServer.R"))
source(file.path("QualimapModule", "qualimapModuleServer.R"))
source(file.path("ComparisonModule", "comparisonModuleServer.R"))

server <- function(input, output, session) {
  vcfModuleServer("vcf")
  qualimapModuleServer("qualimap")
  comparisonModuleServer("comparison")
}
