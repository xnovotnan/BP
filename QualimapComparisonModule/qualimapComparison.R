library(tidyverse)
library(magrittr)
library(fmsb)
library(scales)
library(patchwork)
options(scipen = 999)

# Process QUALIMAP files
process_qualimap_folders <- function(comparison_folder) {
  results_df <- tibble(
    `sample` = character(),
    `number_of_bases` = numeric(),
    `number_of_contigs` = numeric(),
    `number_of_reads` = numeric(),
    `number_of_mapped_reads` = numeric(),
    `number_of_mapped_paired_reads_both_in_pair` = numeric(),
    `number_of_mapped_paired_reads_singletons` = numeric(),
    `number_of_mapped_bases` = numeric(),
    `number_of_sequenced_bases` = numeric(),
    `number_of_duplicated_reads_flagged` = numeric(),
    `mean_mapping_quality` = numeric(),
    `mean_insert_size` = numeric(),
    `std_insert_size` = numeric(),
    `median_insert_size` = numeric(),
    `mean_coverageData` = numeric(),
    `std_coverageData` = numeric(),
    `number_of_A's` = numeric(),
    `number_of_C's` = numeric(),
    `number_of_T's` = numeric(),
    `number_of_G's` = numeric(),
    `number_of_N's` = numeric(),
    `GC_percentage` = numeric()
  )
  required_columns <- colnames(results_df)[-1]
  
  subdirs <- list.dirs(comparison_folder, full.names = TRUE, recursive = FALSE)
  for (folder in subdirs) {
    incProgress(0.1, detail = "Processing file...")
    file_path <- file.path(folder, "genome_results.txt")
    if (file.exists(file_path)) {
      folder_results <- tibble(sample = basename(folder))
      
      lines <- readLines(file_path)
      lines <- lines[!grepl('^>>>>>', lines)]
      start_index <- grep("number of bases", lines) 
      end_index <- grep("std coverageData", lines)
      lines <- lines[start_index:end_index]
      for (line in lines) {
        if (line == "") next
        parts <- strsplit(line, "=")[[1]]
        key <- trimws(parts[1]) 
        key <- gsub(" ", "_", key)
        key <- gsub("\\(", "", key)
        key <- gsub("\\)", "", key)
        value <- trimws(parts[2]) 
        value <- sub("([0-9,\\.]+).*", "\\1", value)
        value <- gsub(",", "", value) %>%
          as.numeric()
        if (key %in% required_columns) {
          folder_results[[key]] <- value
        }
      }
      results_df <- bind_rows(results_df, folder_results)
    } else {
      stop(paste("File genome_results.txt not found in:", basename(folder)))
    }
  }
  results_df
}

# Reference
bases_contigs_comparison <- function(samples_df, attribute, plot_title) {
  samples_df %<>% mutate(label = paste0(sample, " - \n", label_comma()(get(attribute))))
  
  p <- ggplot(samples_df, aes(x = get(attribute), y = reorder(sample, get(attribute)), fill = sample)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = plot_title) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    geom_text(aes(x = get(attribute) / 2, label = label), 
              size = 4, 
              color = "black") +
    scale_fill_brewer(palette = "Set3")
  p
}

# Read Statistics
reads_comparison <- function(samples_df){
  samples_df %<>% select("sample",
                         "number_of_reads",
                         "number_of_mapped_reads") %>%
    pivot_longer(cols = c("number_of_reads","number_of_mapped_reads"), 
                 names_to = "reads", 
                 values_to = "count")
  
  p <- ggplot(samples_df, aes(x=count, y=reorder(sample, count), fill=reads, color=sample)) +
    geom_col(position = position_dodge(0.8), width = 0.7, linewidth=1) +
    labs(title = "Number of Reads and Mapped Reads",
         x = "Number of Reads",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "right",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=12))+
    scale_color_brewer(palette = "Set3")+ 
    scale_fill_manual(values = c("number_of_reads"="lightblue3", "number_of_mapped_reads"="lightblue1"), 
                      labels=c("number_of_reads"="Total Reads", "number_of_mapped_reads"="Mapped Reads"))+
    scale_x_continuous(labels = comma)+
    guides(color = guide_legend(override.aes = list(fill = "white")))
  p
}
mapped_paired_reads <- function(samples_df){
  p <- ggplot(samples_df, aes(x=number_of_mapped_paired_reads_both_in_pair, 
                              y=reorder(sample, number_of_mapped_paired_reads_both_in_pair), 
                              group = 1))+
    geom_line(color = "green3", linewidth = 1) +
    geom_point(color = "green3", size = 5)+
    labs(title = "Number of Mapped Paired Reads (Both in pair)",
         x = "Number of Reads",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))+
    scale_x_continuous(labels = comma)
  p
}
mapped_paired_reads_singletons <- function(samples_df){
  p <- ggplot(samples_df, aes(x=number_of_mapped_paired_reads_singletons, 
                              y=reorder(sample, number_of_mapped_paired_reads_singletons), 
                              group = 1))+
    geom_line(color = "purple2", linewidth = 1) +
    geom_point(color = "purple2", size = 5)+
    labs(title = "Number of Mapped Paired Reads (Singletons)",
         x = "Number of Reads",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))+
    scale_x_continuous(labels = comma)
  p
}
mapped_bases_comparison <- function(samples_df){
  samples_df %<>% select("sample",
                         "number_of_sequenced_bases",
                         "number_of_mapped_bases") %>%
    pivot_longer(cols = c("number_of_sequenced_bases",
                          "number_of_mapped_bases"), 
                 names_to = "reads", values_to = "count")

  p <- ggplot(samples_df, aes(x=count, y=reorder(sample, count), fill=reads, color=sample)) +
    geom_col(position = position_dodge(0.8), width = 0.7, linewidth=1) +
    labs(title = "Number of Mapped and Sequenced Bases (bp)",
         x = "Number of Bases",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "right",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=12))+
    scale_color_brewer(palette = "Set3")+ 
    scale_fill_manual(values = c("number_of_mapped_bases"="lightblue3", "number_of_sequenced_bases"="lightblue1"), 
                      labels=c("number_of_mapped_bases"="Total Bases", "number_of_sequenced_bases"="Mapped Bases"))+
    scale_x_continuous(labels = comma)+
    guides(color = guide_legend(override.aes = list(fill = "white")))
  p
}
duplicated_reads <- function(samples_df){
  p <- ggplot(samples_df, aes(x=number_of_duplicated_reads_flagged, 
                              y=reorder(sample, number_of_duplicated_reads_flagged), 
                              group = 1))+
    geom_line(color = "red3", linewidth = 1) +
    geom_point(color = "red3", size = 5)+
    labs(title = "Number of Duplicated Reads",
         x = "Number of Reads",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))+
    scale_x_continuous(labels = comma)
  p
}
mapping_quality_comparison <- function(samples_df){
  p <- ggplot(samples_df, aes(x=mean_mapping_quality, 
                              y=reorder(sample, mean_mapping_quality), 
                              color=sample)) +
    geom_boxplot(width=0.5) +
    labs(title = "Mean Mapping Quality",
         x = "Mean Quality",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))+
    scale_fill_brewer(palette = "Set3")
  p
}

# Insert Size
insert_size_comparison <- function(samples_df){
  samples_df %<>% mutate(label = paste0(sample, " - \n", label_comma()(mean_insert_size)))
  p <- ggplot(samples_df, aes(x = mean_insert_size, y = reorder(sample, mean_insert_size) , fill = sample, group = 1)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 5, show.legend = FALSE)+
    labs(title = "Histogram of Insert Size", x="Insert Size", y="Sample") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    geom_text(aes(x = mean_insert_size / 2, label = label), 
              size = 4, 
              color = "black") +
    scale_fill_brewer(palette = "Set3")
  p  
}
insert_size_mean_comparison <- function(samples_df){
  p <- ggplot(samples_df, aes(y = reorder(sample, mean_insert_size))) +
    geom_errorbar(aes(xmin = max(mean_insert_size - std_insert_size, 0),
                      xmax = mean_insert_size + std_insert_size), 
                  width = 0.2, color = "black") +  
    geom_point(aes(x = mean_insert_size, color = "Mean", shape = "Mean"), size = 5) +  
    geom_point(aes(x = median_insert_size, color = "Median", shape = "Median"), size = 5) + 
    theme_minimal() +
    labs(title = "Distribution of Insert Size",
         x = "Insert Size", y = element_blank(), 
         color = "Measure", shape = "Measure") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
          legend.position = "bottom", legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=12)) +
    scale_color_manual(values = c("Mean" = "blue", "Median" = "red")) + 
    scale_shape_manual(values = c("Mean" = 16, "Median" = 17)) +
    scale_x_continuous(labels = comma)
  p
}

# Data Coverage
coverage_comparison <- function(samples_df){
  samples_df %<>% mutate(label = paste0(sample, " - \n", label_comma()(mean_coverageData)))
  p <- ggplot(samples_df, aes(x = mean_coverageData, y = reorder(sample, mean_coverageData) , fill = sample, group = 1)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_line(color = "black", linewidth = 1) +
    geom_point(color = "black", size = 5, show.legend = FALSE)+
    labs(title = "Histogram of Data Coverage", x="Data Coverage", y="Sample") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    geom_text(aes(x = mean_coverageData / 2, label = label), 
              size = 4, 
              color = "black") +
    scale_fill_brewer(palette = "Set3")
 p
}
coverage_mean_comparison <- function(samples_df){
  p <- ggplot(samples_df, aes(x=mean_coverageData, y=reorder(sample, mean_coverageData))) +
    geom_errorbar(aes(xmin = mean_coverageData - std_coverageData, 
                      xmax = mean_coverageData + std_coverageData), 
                  width = 0.2, color = "red") + 
    geom_point(size=5, color="orange") +
    labs(title = "Distribution of Data Coverage",
         x = "Data Coverage",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}

# ACTG Content
gc_comparison <- function(samples_df){
  p <- ggplot(samples_df, aes(x=GC_percentage, 
                              y=reorder(sample, GC_percentage), 
                              group = 1))+
    geom_line(color = "blue2", linewidth = 1) +
    geom_point(color = "blue2", size = 5)+
    labs(title = "GC Percentage (%)",
         x = "GC Content",
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}
stacked_actg <- function(samples_df){
  samples_df %<>% select(
    "sample",
    "A" = "number_of_A's",
    "C" = "number_of_C's",
    "T" = "number_of_T's",
    "G" = "number_of_G's",
    "N" = "number_of_N's"
  )
  
  samples_df %<>%
    pivot_longer(cols = -sample, names_to = "base", values_to = "count") %>%
    group_by(sample) %>%
    mutate(percentage = round(count *100/ sum(count), 2))
  samples_df$base <- with(samples_df, reorder(base, percentage, FUN = sum))
  
  p <- ggplot(samples_df, aes(x = percentage, y = sample, fill = base)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage",
         y = element_blank(), 
         title= "ACTG Distribution (%)") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "right",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size=12))
  p
}

