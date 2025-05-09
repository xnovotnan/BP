library(shiny)
source(file.path("VCFModule", "vcfModuleServer.R"))
source(file.path("VCFComparisonModule", "vcfComparisonServer.R"))
source(file.path("QualimapModule", "qualimapModuleServer.R"))
source(file.path("QualimapComparisonModule", "qualimapComparisonServer.R"))
source(file.path("FastqcModule", "fastqcModuleServer.R"))

# The server function defines the main server logic in the application.
# It connects all modules via their IDs. 

server <- function(input, output, session) {
  vcfModuleServer("vcf")
  qualimapModuleServer("qualimap")
  comparisonModuleServer("qualimapComparison")
  vcfComparisonServer("vcfComparison")
  fastqcModuleServer("fastqc")
  
  # deletes temporary files
  shiny::onStop(function() {
    unlink("Temp/*", recursive = TRUE)
  })
}
