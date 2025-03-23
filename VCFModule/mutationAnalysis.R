library(tidyverse)
library(magrittr)
library(patchwork)
library(plotly)
library(ggridges)
library(scales)
options(scipen = 999)

# -----------------------------------
# TYPY MUTACII A ICH DISTRIBUCIA
# -----------------------------------

# Mutation Counts
mutation_donut <- function(vcfTibble, subtypes=FALSE, valueType){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  mutationCounts <- vcfTibble %>%
    group_by(TYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(count * 100 / sum(count), 2),
           ymax = cumsum(percentage),
           ymin = c(0, head(ymax, n=-1)),
           labelPosition = (ymax + ymin) / 2,
           label = case_when(valueType == "Percentage" ~ paste0(TYPE, "s\n", percentage, "%"), 
                             TRUE ~ paste0(TYPE, "s\n", label_comma()(count))))
           
  p <- ggplot(mutationCounts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=TYPE)) +
    geom_rect(show.legend = FALSE) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    scale_fill_brewer(palette = "Set3") +
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=3) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    labs(title = "Variant Distribution")
  p
}

mutation_distribution <- function(vcfTibble, subtypes=FALSE, valueType){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  mutationCounts <- vcfTibble %>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    group_by(CHROM) %>%
    mutate(percentage = count * 100/ sum(count))
  if (valueType == "Percentage") {mutationCounts$count <- mutationCounts$percentage}

  p <- ggplot(mutationCounts, aes(x = CHROM, y = count, fill = TYPE)) +
    geom_bar(stat = "identity", position = "stack", show.legend = FALSE) +
    labs(x = NULL,y = NULL, title = "Variant Distribution on Chromosomes") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")
  p
}


# SNV Analysis
snv_types <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV")
  mutationCounts <- vcfTibble %>%
    group_by(SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = count / sum(count),
           ymax = cumsum(percentage),
           ymin = c(0, head(ymax, n=-1)),
           labelPosition = (ymax + ymin) / 2,
           label =paste0(SUBTYPE, "s\n", label_comma()(count),"\n", round(percentage * 100,2), "%"))
  
  p <- ggplot(mutationCounts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=SUBTYPE)) +
    geom_rect(show.legend = FALSE) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=3) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
    labs(title = "SNV Types Distribution") +
    scale_fill_brewer(palette = "Set3")
  p
}

snv_class_barplot <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV") %>% 
    count(SNV_TYPE = paste(REF, ALT, sep = ">")) %>% 
    mutate(percentage = round(n *100/ sum(n), 2),
           label = paste0(SNV_TYPE, "\n", percentage, "%"))
  vcfTibble$SNV_TYPE <- with(vcfTibble, reorder(SNV_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcfTibble, aes(x = percentage, y = reorder(SNV_TYPE, n), fill = SNV_TYPE)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "SNV Types Percentages", x = NULL, y = "Percentage") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    geom_text(aes(x = percentage/2, label = label), size = 3, color = "black")+
    scale_fill_brewer(palette = "Set3")
  
  p
}

snv_class_boxplot <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV") %>% group_by(CHROM) %>%
    count(SNV_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2),
           subtype = ifelse(SNV_TYPE %in% c("A>G", "G>A", "C>T", "T>C"), 
                            "Transition", "Transversion"))
  vcfTibble$SNV_TYPE <- with(vcfTibble, reorder(SNV_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcfTibble, aes(x = reorder(SNV_TYPE, -percentage), y = percentage, 
                             color = subtype, fill = SNV_TYPE)) +
    geom_boxplot() +
    labs(title = "Distribution of SNV Type",
         x = element_blank(),
         y = "Percentage") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.title = element_blank()) +
    scale_fill_brewer(palette = "Set3")
  
  p
}

snv_class_stacked <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV") %>% group_by(CHROM) %>%
    count(SNV_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  vcfTibble$SNV_TYPE <- with(vcfTibble, reorder(SNV_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcfTibble, aes(x = percentage, y = CHROM, fill = SNV_TYPE)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage", y = element_blank(), title= "SNV Types Distribution on Chromosomes") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "right",  plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.title = element_blank())
  
  p
}

snv_values <- function(vcfTibble){
  values <- vcfTibble %>% group_by(TYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2)) %>%
    filter(TYPE == "SNV")
}


# INDEL Analysis
indel_values <- function(vcfTibble){
  values <- vcfTibble %>% group_by(TYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2)) %>%
    filter(TYPE == "INDEL")
}

indel_types <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL")
  mutationCounts <- vcfTibble %>%
    group_by(SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = count / sum(count),
           ymax = cumsum(percentage),
           ymin = c(0, head(ymax, n=-1)),
           labelPosition = (ymax + ymin) / 2,
           label =paste0(SUBTYPE, "s\n", label_comma()(count),"\n" ,round(percentage * 100,2), "%"))
  
  p <- ggplot(mutationCounts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=SUBTYPE)) +
    geom_rect(show.legend = FALSE) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=3) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
    labs(title = "Insertion-Deletion Distribution") +
    scale_fill_brewer(palette = "Set3")
  p
}

indel_stacked <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>% group_by(CHROM, SUBTYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2))

  p <- ggplot(vcfTibble, aes(x = percentage, y = CHROM, fill = SUBTYPE)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage", y = element_blank(), title= "Insertions and Deletions on Chromosomes") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "right",  plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          legend.title = element_blank())
  p
}

indel_length <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  
  p <- ggplot(vcfTibble, aes(x = length)) +
    geom_histogram(binwidth = 1, position = "dodge") +
    labs(x = "INDEL length (bp)", y = "Mutation Count") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")
  p
}

indel_length_boxplot <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF))) %>%
    group_by(CHROM, SUBTYPE)
 
  p <- ggplot(vcfTibble, aes(x = SUBTYPE, y = length, fill = SUBTYPE)) +
    geom_boxplot(show.legend = FALSE) +
    labs(x = element_blank(), y = "Length (bp)") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")
  p
}

indel_length_avg <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  mean(vcfTibble$length)
}

indel_length_med <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  median(vcfTibble$length)
}



