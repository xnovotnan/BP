library(tidyverse)
library(magrittr)
library(patchwork)
library(plotly)
library(ggridges)
library(scales)
theme_set(theme_minimal() + 
            theme(
              plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
              legend.position = "bottom"
            ))
options(scipen = 999)
my_colors <- c("skyblue", "salmon", "orange", "lightgreen", "purple", "yellow","lightpink")


# -----------------------------------
# TYPY MUTACII A ICH DISTRIBUCIA
# -----------------------------------

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
    labs(title = "Variant summary")
  p
}

mutation_distribution <- function(vcfTibble, subtypes=FALSE, valueType){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  mutationCounts <- vcfTibble %>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    group_by(CHROM) %>%
    mutate(percentage = count * 100/ sum(count))
  if (valueType == "Percentage") {mutationCounts$count <- percentage}

  p <- ggplot(mutationCounts, aes(x = CHROM, y = count, fill = TYPE)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL,y = NULL, title = "Variant Distribution on Chromosomes") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")
  p
}

snv_types <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV")
  mutationCounts <- vcfTibble %>%
    group_by(SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = count / sum(count),
           ymax = cumsum(percentage),
           ymin = c(0, head(ymax, n=-1)),
           labelPosition = (ymax + ymin) / 2,
           label =paste0(SUBTYPE, "s\n", round(percentage * 100,2), "%"))
  
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

# SNV types
snv_classes <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV") %>% group_by(SUBTYPE) %>%
    count(SNV_TYPE = paste(REF, ALT, sep = ">")) %>% 
    ungroup() %>%
    mutate(percentage = round(n *100/ sum(n), 2),
           label = paste0(SNV_TYPE, "\n", percentage, "%"))
  
  transitions <- vcfTibble %>% filter(SUBTYPE == "Transition")
  transversions <- vcfTibble %>% filter(SUBTYPE == "Transversion")
  
  p1 <- ggplot(transitions, aes(x = percentage, y = reorder(SNV_TYPE, n), fill = SNV_TYPE)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Transitions", x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    geom_text(aes(x = percentage / 2, label = label), size = 3, color = "black") +
    scale_fill_brewer(palette = "Set3")
  
  
  p2 <- ggplot(transversions, aes(x = percentage, y = reorder(SNV_TYPE, n), fill = SNV_TYPE)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Transversions", x = NULL, y = NULL) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    geom_text(aes(x = percentage/2, label = label), size = 3, color = "black")+
    scale_fill_brewer(palette = "Set3")
  
  (p1 | p2)
}

snv_class_combined <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV") %>% group_by(CHROM, SUBTYPE) %>%
    count(SNV_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  
  p1 <- ggplot(vcfTibble, aes(x = reorder(SNV_TYPE, -percentage), y = percentage, fill = SUBTYPE)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = "Distribution of SNV Type",
         x = element_blank(),
         y = "Percentage") +
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) +
    scale_fill_brewer(palette = "Set3")
  
  p2 <- ggplot(vcfTibble, aes(x = SUBTYPE, y = percentage, fill = SUBTYPE)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = NULL,
         x = element_blank(),
         y = element_blank()) +
    theme_minimal()
  (p1 | p2)
}


snv_class_stacked <- function(vcfTibble){
  vcfTibble %<>% filter(TYPE == "SNV") %>% group_by(CHROM, SUBTYPE) %>%
    count(SNV_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  transitions <- vcfTibble %>% filter(SUBTYPE == "Transition")
  transversions <- vcfTibble %>% filter(SUBTYPE == "Transversion")
  
  p1 <- ggplot(transitions, aes(x = CHROM, y = percentage, fill = SNV_TYPE)) +
    geom_bar(stat = "identity", position = "stack", show.legend = FALSE) + 
    labs(x = element_blank(), y = element_blank(), fill = "SNV Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_fill_brewer(palette = "Set3") 
  
  p2 <- ggplot(transversions, aes(x = CHROM, y = percentage, fill = SNV_TYPE)) +
    geom_bar(stat = "identity", position = "stack", show.legend = FALSE) + 
    labs(x = element_blank(), y = element_blank(), fill = "SNV Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    scale_fill_brewer(palette = "Set3") 
  
  (p1 | p2)
}