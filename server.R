library(shiny)
source(file.path("VCFModule", "vcfModuleServer.R"))
source(file.path("VCFComparisonModule", "vcfComparisonServer.R"))
source(file.path("QualimapModule", "qualimapModuleServer.R"))
source(file.path("QualimapComparisonModule", "qualimapComparisonServer.R"))
source(file.path("FastqcModule", "fastqcModuleServer.R"))

server <- function(input, output, session) {
  vcfModuleServer("vcf")
  qualimapModuleServer("qualimap")
  comparisonModuleServer("qualimapComparison")
  vcfComparisonServer("vcfComparison")
  fastqcModuleServer("fastqc")
  
  shiny::onStop(function() {
    unlink("Temp/*", recursive = TRUE)
  })
}


