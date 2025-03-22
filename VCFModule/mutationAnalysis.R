library(tidyverse)
library(magrittr)
library(patchwork)
library(plotly)
library(ggridges)
theme_set(theme_minimal() + 
            theme(plot.title = element_text(hjust = 0.5, size = 10, face = "bold")))
options(scipen = 999)
my_colors <- c("skyblue", "salmon", "orange", "lightgreen", "purple")

# -----------------------------------
# TYPY MUTACII A ICH DISTRIBUCIA
# -----------------------------------

# Distribúcia mutácií - barplot a pie chart
mut_summary <- function(vcfTibble, subtypes=FALSE){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  
  mutationCounts <- vcfTibble %>% 
    group_by(TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    mutate(percentage = count * 100/ sum(count)) 
  
  p <-   ggplot(mutationCounts, aes(x = TYPE, y = percentage, fill = TYPE)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme(legend.position="none") +
    labs(y=NULL, x=NULL, title = "Variant Summary")+
    scale_fill_manual(values = my_colors)
  p
}


# Distribúcia mutácií na jednotlivých chromozómoch 
# možnosť zvoliť úroveň klasifikácie mutácií aj konkrétny chromozóm
mut_dist <- function(vcfTibble, subtypes=FALSE, valueType){
  if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
  
  mutationCounts <- vcfTibble %>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    group_by(CHROM) %>%
    mutate(percentage = count * 100/ sum(count))

  
  if (valueType == "Percentage") {
    p <- ggplot(mutationCounts, aes(x = CHROM, y = percentage, fill = TYPE)) +
      geom_bar(stat = "identity", position = "stack")
  }
  else{
    p <- ggplot(mutationCounts, aes(x = CHROM, y = count, fill = TYPE)) +
      geom_bar(stat = "identity", position = "stack")
  }

  p <- p +
    labs(x = NULL,
         y = NULL,
         title = "Variant Distribution on Chromosomes") +
    coord_flip() +
    scale_fill_manual(values = my_colors)
  
  p
}


# # Typy mutácií na kruhovom grafe 
# mutation_summary_circle <- function(vcfTibble, subtypes=FALSE){
#   if (subtypes) vcfTibble$TYPE <- vcfTibble$SUBTYPE
#   mutation_summary <- vcfTibble %>% 
#     group_by(CHROM, TYPE) %>% 
#     summarise(count = n(), .groups = "drop")
#   
#   ggplot(mutation_summary, aes(x = CHROM, y = count, fill = TYPE)) +
#     geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
#     theme_minimal() +
#     labs(title = "Mutation types", 
#          y = "Number of mutations")+
#     coord_polar()+
#     scale_y_log10()
# }
# 
# 
# # Distribúcia SNP typov
# SNV_types <- function(vcfTibble){
#   vcfTibble %<>% filter(TYPE == "SNP") %>%  
#     count(SNP_TYPE = paste(REF, ALT, sep = ">"))
#   
#   p <- ggplot(vcfTibble, aes(x = reorder(SNP_TYPE, -n), y = n, fill = SNP_TYPE)) +
#     geom_bar(stat = "identity", color = "black") +
#     theme_minimal() +
#     labs(title = "SNP Mutations Distribution", y=NULL, x=NULL) +
#     theme(legend.position = "none")+
#     coord_flip()
#   
#   ggplotly(p)
# }