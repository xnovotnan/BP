library(VariantAnnotation)
library(BiocManager)
library(tidyverse)
library(magrittr)

setwd("/Users/macbook/Documents/BP/grafy/qualimap")
path <- "/Users/macbook/Documents/BP/data/qualimap_06"
files <- list.files(path = path, full.names = TRUE)

for (file in files){
  data <- read.delim(file, header = TRUE, sep = "\t")
  if (grepl("mapped_reads_nucleotide_content", file) | grepl("homopolymer_indels", file)){
    print(file)
  }
  else{
    file_name <- basename(file)
    file_name <- substr(file_name, 1, nchar(file_name) - 4)

    create_plots(data, file_name)
  }
}


create_plots <- function(data, file_name){
  columns <- colnames(data)
  
  histogram <- ggplot(data, aes_string(x = columns[1], y = columns[2])) +
    geom_bar(stat = "identity", fill = "blue", color = "black") +
    labs(
      title = paste(file_name),
      x = columns[1],
      y = columns[2]
    ) +
    theme_minimal()
  
  ggsave(paste(file_name, "_histogram.png"), plot = histogram, width = 8, height = 6, dpi = 300)
  
  linechart <- ggplot(data, aes_string(x = columns[1], y = columns[2])) +
    geom_line(color="blue") +
    labs(
      title = paste(file_name),
      x = columns[1],
      y = columns[2]
    ) +
    theme_minimal()

  ggsave(paste(file_name, "_linechart.png"), plot = linechart, width = 8, height = 6, dpi = 300)
  print("SAVED")
}

