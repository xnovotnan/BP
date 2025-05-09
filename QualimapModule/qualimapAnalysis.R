library(tidyverse)
library(magrittr)
options(scipen = 999)

# The get_genome_results function reads and returns the content of the genome_results.txt file 
# from the selected QUALIMAP folder, stopping with an error if the file is missing.
get_genome_results <- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "genome_results.txt")
  if (!file.exists(file_path)) {
    stop("File genome_results.txt not found in this folder.")
  }
  readLines(file_path)
}

# The process_qualimap function extract relevant data from the QUALIMAP output 
# and creates a dictionary with these values. 
process_qualimap <- function(qualimap_folder) {
  incProgress(0.1, detail = "Reading file...")
  lines <- get_genome_results(qualimap_folder)
  lines <- lines[!grepl('^>>>>>', lines)]
  start_index <- grep("bam file", lines) 
  end_index <- grep("std coverageData", lines)
  lines <- lines[start_index:end_index]
  results <- list()
  incProgress(0.5, detail = "Processing data...")
  for (line in lines) {
    if (line == "") next
    parts <- strsplit(line, "=")[[1]]
    key <- trimws(parts[1])
    value <- trimws(parts[2])
    results[[key]] <- value
  }
  results
}

# Data Coverage
process_qualimap_coverage <- function(qualimap_folder){
  # reads data related to coverage and processes it
  lines <- get_genome_results(qualimap_folder)
  start_index <- grep("std coverageData", lines) + 1
  end_index <- grep("Coverage per contig", lines) - 1
  lines <- lines[start_index:end_index]
  percentages <- numeric()
  
  for (line in lines) {
    if (line == "") next
    line <- trimws(line)
    parts <- strsplit(line, " ")[[1]]
    percentage <- parts[[4]] 
    percentage <- substr(percentage, 1, nchar(percentage) - 1) %>%  as.numeric(percentage)
    percentages <- c(percentages, percentage)
  }
  df <- tibble(X = seq_along(percentages), Y = percentages)
  
  ggplot(df, aes(x = X, y = Y)) +
    geom_bar(stat = "identity", fill = "skyblue", alpha = 0.5) +
    geom_line(color = "blue", size = 1) + 
    geom_point(color = "blue", size = 2) + 
    labs(x = "Coverage (X)", y = "Percentage (%)", title = "Genome Fraction Coverage")+
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
      )
}
process_qualimap_coverage_pc <- function(qualimap_folder) {
  # reads data related to coverage per chromosome and processes it
  lines <- get_genome_results(qualimap_folder)
  start_index <- grep("chr1\t", lines)
  end_index <- grep("chrY", lines)
  lines <- lines[start_index:end_index]
  
  df <- read.table(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chromosome", "contig_size", "mapped_reads", "mean_coverage", "std_coverage")
  df$chromosome <- factor(df$chromosome, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
  
  p <- ggplot(df, aes(x = chromosome, y = mean_coverage, group = 1)) +
    geom_line(color = "salmon", size = 1) +
    geom_point(color = "salmon", size = 2) +
    labs(title = "Mean Coverage Across Chromosomes", x = element_blank(), y = "Mean coverage")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12))
  p
}

# ACTG
extract_value <- function(value) {
  # extracts data related to ACTG content 
  value <- sub("([0-9,]+).*", "\\1", value)  
  value <- gsub(",", "", value) %>%
    as.numeric()
  return(value)
}
process_ACTG_content <- function(data){
  # reads data related to ACTG content and processes it
  data <- data[c("number of A's", "number of C's", "number of T's", "number of G's", "number of N's")]
  values <- sapply(data, extract_value)
  labels <- c("A", "C", "T", "G", "N")
  df <- data.frame(Base = labels, Count = values)
  df$Percentage <- df$Count / sum(df$Count) * 100
  
  p <- ggplot(df, aes(x = reorder(Base, -Percentage), y = Percentage)) +
    geom_bar(stat = "identity", fill = "lightgreen", alpha=0.5) +
    labs(title = "Base Pair Counts", 
         x = element_blank(), 
         y = "Percentage (%)")+
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  p
}

# The find_qualimap_png function locates the .png file with name from the parameters   
# in the selected QUALIMAP folder, returnign NULL if the file is missing.
find_qualimap_png <- function(qualimap_folder, name){
  file_path <- file.path(qualimap_folder, paste("images_qualimapReport", name, sep="/"))
  if (!file.exists(file_path)) {
    return(NULL)
  }
  file_path
}

