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
mutation_donut <- function(vcfTibble, subtypes=FALSE){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  mutationCounts <- vcfTibble %>%
    group_by(TYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(count * 100 / sum(count), 2),
           ymax = cumsum(percentage),
           ymin = c(0, head(ymax, n=-1)),
           labelPosition = (ymax + ymin) / 2,
           label = paste0(TYPE, "s\n", percentage, "%\n", label_comma()(count)))
  
  p <- ggplot(mutationCounts, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=TYPE)) +
    geom_rect(show.legend = FALSE) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    scale_fill_brewer(palette = "Set3") +
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=4) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    labs(title = "Variant Distribution")
  p
}
mutation_distribution <- function(vcfTibble, subtypes=FALSE){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  mutationCounts <- vcfTibble %>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    group_by(CHROM) %>%
    mutate(percentage = count * 100/ sum(count))
  
  p <- ggplot(mutationCounts, aes(x = CHROM, y = percentage, fill = TYPE)) +
    geom_bar(stat = "identity", position = "stack", show.legend = FALSE) +
    labs(x = NULL, y = NULL, title = "Variant Distribution on Chromosomes (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}

# SNP Analysis
snp_value <- function(vcfTibble){
  values <- vcfTibble %>% group_by(TYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2)) %>%
    filter(TYPE == "SNP")
}
snp_types_donut <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNP")
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
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=4) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
    labs(title = "SNP Types Distribution") +
    scale_fill_brewer(palette = "Set3")
  p
}
snp_class_stacked <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNP") %>% group_by(CHROM) %>%
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  vcfTibble$SNP_TYPE <- with(vcfTibble, reorder(SNP_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcfTibble, aes(x = percentage, y = CHROM, fill = SNP_TYPE)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage", 
         y = element_blank(), 
         title= "SNP Types Distribution on Chromosomes") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "right",)
  p
}
snp_class_boxplot <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNP") %>% group_by(CHROM) %>%
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2),
           subtype = ifelse(SNP_TYPE %in% c("A>G", "G>A", "C>T", "T>C"), 
                            "Transition", "Transversion"))
  vcfTibble$SNP_TYPE <- with(vcfTibble, reorder(SNP_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcfTibble, 
              aes(x = reorder(SNP_TYPE, -percentage), 
                  y = percentage, 
                  color = subtype, 
                  fill = SNP_TYPE)) +
    geom_boxplot() +
    labs(title = "Distribution of SNP Type",
         x = element_blank(),
         y = element_blank()) +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    scale_fill_brewer(palette = "Set3")
  p
}
snp_class_barplot <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNP") %>% 
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>% 
    mutate(percentage = round(n *100/ sum(n), 2),
           label = paste0(SNP_TYPE, "\n", percentage, "%"))
  vcfTibble$SNP_TYPE <- with(vcfTibble, reorder(SNP_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcfTibble, aes(x = percentage, y = reorder(SNP_TYPE, n), fill = SNP_TYPE)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "SNP Types Percentages (%)", 
         x = NULL, 
         y = "Percentage") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
    geom_text(aes(x = percentage/2, label = label), size = 4, color = "black")+
    scale_fill_brewer(palette = "Set3")
  p
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
    geom_text(x=3.5, aes(y=labelPosition, label=label), size=4) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))+
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
    labs(x = element_blank(), 
         y = element_blank(), 
         title= "Insertions and Deletions on Chromosomes (%)") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank())
  p
}
indel_length_avg <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  round(mean(vcfTibble$length),2)
}
indel_length_med <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  median(vcfTibble$length)
}
indel_length <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  
  p <- ggplot(vcfTibble, aes(x = length)) +
    geom_histogram(binwidth = 1, position = "dodge", color="purple") +
    labs(x = "INDEL length (bp)", 
         y = "Mutation Count",
         title = "INDEL Length") +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}
indel_length_boxplot <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF))) %>%
    group_by(CHROM, SUBTYPE)
  
  p <- ggplot(vcfTibble, aes(x = SUBTYPE, y = length, fill = SUBTYPE)) +
    geom_boxplot(show.legend = FALSE) +
    labs(x = element_blank(), 
         y = "Length (bp)") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3")+
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}

# Quality
quality_avg <- function(vcfTibble){
  round(mean(vcfTibble$QUAL),2)
}
quality_med <- function(vcfTibble){
  median(vcfTibble$QUAL)
}
quality_bar <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x=QUAL)) +
    geom_bar(color="lightgreen",alpha=0.5) +
    labs(title = "Quality across the file", 
         x = "Position", 
         y="Quality") +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  p
}
quality_on_chroms <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x=CHROM, y=QUAL)) +
    geom_violin(color="orange") +
    labs(title = "Quality across the chromosomes",
         x = element_blank(),
         y = "Quality") +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size=12),
          axis.text.y = element_text(size = 12))
  p
}

# Allele Frequency
allele_freq_avg <- function(vcfTibble){
  round(mean(vcfTibble$AF),2)
}
allele_freq_med <- function(vcfTibble){
  median(vcfTibble$AF)
}
allele_freq_hexbin <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x=seq_along(AF), y = AF)) +
    geom_hex() +
    labs(title = "Allele Frequency across the file", 
         x = "Position", 
         y="Allele Frequency") +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  p
}
allele_freq_on_chroms <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x=AF, y=CHROM)) +
    geom_density_ridges(fill="lightpink") +
    labs(title = "Allele Frequency across the chromosomes",
         x= element_blank(),
         y=element_blank()) +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}

# Read Depth
read_depth_avg <- function(vcfTibble){
  round(mean(vcfTibble$DP),2)
}
read_depth_med <- function(vcfTibble){
  median(vcfTibble$DP)
}
read_depth_density <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x=DP)) +
    geom_bar(color="lightgreen",alpha=0.5) +
    labs(title = "Read Depth across the file", 
         x = "Position", 
         y="Read Depth") +
    theme_minimal() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}
read_depth_on_chroms <- function(vcfTibble){
  p <- ggplot(vcfTibble, aes(x=CHROM, y=DP)) +
    geom_violin(color="orange") +
    labs(title = "Read depth across chromosomes",
         x = "Chromosomes",
         y = "Read depth") +
    theme_minimal() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size=12),
          axis.text.y = element_text(size = 12))
  p
}







