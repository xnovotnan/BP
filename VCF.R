library(tidyverse)
library(magrittr)
library(patchwork)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(plotly)
library(circlize)
library(ComplexHeatmap)

setwd("/Users/macbook/Documents/BP")
vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N.vcf"



# -----------------------------------
# GRAFY KVALITY A HLBKY CITANI
# -----------------------------------

# Hexbin QUAL vs Read depth
qual_vs_rd <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x = DP, y = QUAL)) +
    geom_hex() + 
    #scale_x_log10() +
    theme_minimal()
  ggplotly(p)
}
qual_vs_rd(vcf_complete)


# -----------------------------------
# GRAFY ALELICKEJ FREKVENCIE
# -----------------------------------

# Graf alelickej frekvencie - density plot 
plot_frequency <- function(vcfTibble) {
  p <- ggplot(vcfTibble, aes(x=AF)) +
    geom_density(fill = "skyblue", alpha=0.5) +
    labs(title = "Allele Frequency",
         x = "Mutation",
         y = "Frequency") +
    theme_minimal()
  ggplotly(p)
}
plot_frequency(vcf_complete)


