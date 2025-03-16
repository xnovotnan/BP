library(tidyverse)
library(magrittr)
library(patchwork)
library(plotly)
library(ggridges)
library(ComplexHeatmap)
theme_set(theme_minimal() + 
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.position = "bottom",
              axis.text.x = element_text(angle = 90, hjust = 1)
            ))
options(scipen = 999)
my_colors <- c("skyblue", "salmon", "orange", "lightgreen", "purple")

# -----------------------------------
# PREDSPRACOVANIE DAT, ZAKLADNA ANALYZA A VIZUALIZACIA QUALIMAP SUBORU
# -----------------------------------

# Funkcia get_genome_results najde v qualimap foldri genome_results.txt
get_genome_results <- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "genome_results.txt")
  if (!file.exists(file_path)) {
    stop("File genome_results.txt does not found in this folder.")
  }
  readLines(file_path)
}

# Funkcia process_qualimap() sluzi na spracovanie qualimap suboru - vyberie
# riadky s dôležitými analýzami a spracuje ich do dictionary
process_qualimap <- function(qualimap_folder) {
  lines <- get_genome_results(qualimap_folder)
  lines <- lines[!grepl('^>>>>>', lines)]
  start_index <- grep("number of bases", lines) 
  end_index <- grep("std coverageData", lines)
  lines <- lines[start_index:end_index]
  results <- list()
  
  for (line in lines) {
    if (line == "") next
    parts <- strsplit(line, "=")[[1]]
    key <- trimws(parts[1])
    value <- trimws(parts[2])
    results[[key]] <- value
  }
  results
}


# Funkcia process_qualimap_coverage sluzi na spracovanie coverage dat do grafu
process_qualimap_coverage <- function(qualimap_folder){
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
  p <- ggplot(df, aes(x = X, y = Y)) +
    geom_bar(stat = "identity", fill = "skyblue", alpha = 0.5) +
    geom_line(color = "blue", size = 1) + 
    geom_point(color = "blue", size = 2) + 
    labs(x = "Coverage (X)", y = "Percentage (%)", title = "Coverage Linechart")
  p
}

# Funkcia process_qualimap_coverage_pc sluzi na spracovanie coverage per contig dat
process_qualimap_coverage_pc <- function(qualimap_folder) {
  lines <- get_genome_results(qualimap_folder)
  start_index <- grep("chr1\t", lines)
  end_index <- grep("chrY", lines) #vynechany chrM
  lines <- lines[start_index:end_index]
  
  df <- read.table(text = paste(lines, collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("chromosome", "contig_size", "mapped_reads", "mean_coverage", "std_coverage")
  df$chromosome <- factor(df$chromosome, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
  
  p <- ggplot(df, aes(x = chromosome, y = mean_coverage, group = 1)) +
    geom_line(color = "salmon", size = 1) +
    geom_point(color = "salmon", size = 2) +
    labs(title = "Mean Coverage Across Chromosomes", x = "Chromosome", y = "Mean Coverage") 
  p
}


extract_value <- function(value) {
  value <- sub("([0-9,]+).*", "\\1", value)  
  value <- gsub(",", "", value) %>%
    as.numeric()
  return(value)
}

# Funkcia process_ACTG_content spracuje ACTG content do grafu 
process_ACTG_content <- function(data){
  data <- data[c("number of A's", "number of C's", "number of T's", "number of G's", "number of N's")]
  values <- sapply(data, extract_value)
  labels <- c("A", "C", "T", "G", "N")
  df <- data.frame(Base = labels, Count = values)
  df$Percentage <- df$Count / sum(df$Count) * 100
  
  p <- ggplot(df, aes(x = Base, y = Percentage)) +
    geom_bar(stat = "identity", fill = "salmon", alpha=0.5) +
    labs(title = "Base Pair Counts (Percentage)", x = "Bases", y = "Percentage (%)")
    # scale_fill_manual(values = c("A" = "salmon", "C" = "skyblue", "T" = "orange", "G"= "lightgreen", "N" = "purple"))
  p
}

# Funkcia show_gc_content najde cestu k gc content grafu
show_gc_content <- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "images_qualimapReport/genome_gc_content_per_window.png")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  file_path
}

# Funkcia show_insert_size_across_reference najde cestu k insert size across reference 
show_insert_size_across_reference<- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "images_qualimapReport/genome_insert_size_across_reference.png")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  file_path
}

# Funkcia show_insert_size_histogram najde cestu k insert size histogramu
show_insert_size_histogram <- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "images_qualimapReport/genome_insert_size_histogram.png")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  file_path
}

# Funkcia show_duplication_rate_histogram najde cestu k duplication rate histogramu
show_duplication_rate_histogram <- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "images_qualimapReport/genome_uniq_read_starts_histogram.png")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  file_path
}

# Funkcia show_mapping_quality_histogram najde cestu k mapping quality histogramu
show_mapping_quality_histogram <- function(qualimap_folder){
  file_path <- file.path(qualimap_folder, "images_qualimapReport/genome_mapping_quality_histogram.png")
  if (!file.exists(file_path)) {
    return(NULL)
  }
  file_path
}




