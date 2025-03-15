library(tidyverse)
library(magrittr)
library(patchwork)
library(plotly)
library(ggridges)
library(ComplexHeatmap)
theme_set(theme_minimal() + 
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.position = "bottom"
            ))
options(scipen = 999)

my_colors <- c("skyblue", "salmon", "orange", "lightgreen", "purple")

# -----------------------------------
# TYPY MUTACII A ICH DISTRIBUCIA
# -----------------------------------

# Distribúcia mutácií - barplot a pie chart
mut_summary <- function(vcf_tibble, subtypes=FALSE){
  if (subtypes) vcf_tibble$TYPE <- vcf_tibble$SUBTYPE
  mutation_counts <- vcf_tibble %>% 
    group_by(TYPE) %>% 
    summarise(count = n(), .groups = "drop")
  
  p1 <- ggplot(vcf_tibble, aes(x=TYPE, fill = TYPE)) +
    geom_bar(color = "black") + 
    theme_minimal()+
    coord_flip()+
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    labs(y=NULL, x=NULL)+
    scale_fill_manual(values = my_colors)
  
  
  p2 <- ggplot(vcf_tibble, aes(x=TYPE, fill = TYPE)) +
    geom_bar(color = "black") + 
    theme_minimal()+
    coord_polar()+
    theme(legend.position="none", 
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    labs(y=NULL,x=NULL)+
    scale_fill_manual(values = my_colors)
  
  (p1 + p2) + 
    plot_layout(ncol = 2) + 
    plot_annotation(title = "Variant Summary") &
    theme(plot.title = element_text(hjust = 0.5))
}


# Distribúcia mutácií na jednotlivých chromozómoch 
# možnosť zvoliť úroveň klasifikácie mutácií aj konkrétny chromozóm
mut_dist <- function(vcf_tibble, subtypes=FALSE){
  if (subtypes) vcf_tibble$TYPE <- vcf_tibble$SUBTYPE
  
  p <- ggplot(vcf_tibble, 
              aes(x = CHROM, fill = TYPE)) +
    geom_bar(position = "stack", color = "black") +
    labs(x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
    coord_flip() +
    scale_fill_manual(values = my_colors)
  
  p +
    plot_annotation(title = "Variant Distribution on Chromosomes") &
    theme(plot.title = element_text(hjust = 0.5)) 
}


# Typy mutácií na kruhovom grafe 
mutation_summary_circle <- function(vcf_tibble, subtypes=FALSE){
  if (subtypes) vcf_tibble$TYPE <- vcf_tibble$SUBTYPE
  mutation_summary <- vcf_tibble %>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop")
  
  ggplot(mutation_summary, aes(x = CHROM, y = count, fill = TYPE)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Mutation types", 
         y = "Number of mutations")+
    coord_polar()+
    scale_y_log10()
}


# Distribúcia typov mutácií na lollipop grafe
lollipop <- function(vcf_tibble){ #snp x indel
  mutation_summary <- vcf_tibble %>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    pivot_wider(names_from = TYPE, values_from = count, values_fill = list(count = 0))
  
  colors <- c("SNP" = "green", "INDEL" = "red")
  p <- ggplot(mutation_summary) +
    geom_segment(aes(x=CHROM, xend=CHROM, y=INDEL, yend=SNP), color="grey") +
    geom_point(aes(x=CHROM, y=INDEL,color="INDEL"), size=3, alpha=0.75) +
    geom_point(aes(x=CHROM, y=SNP, color="SNP"), size=3, alpha=0.75) +
    coord_flip()+
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    labs(title = "Number of SNPs and INDELs",
         color = "Legend") +
    scale_color_manual(values = colors)
  ggplotly(p)
}


# Distribúcia subtypov mutácií na lollipop grafe
lollipop_subtypes <- function(vcf_tibble){
  mutation_summary <- vcf_tibble %>% 
    group_by(CHROM, SUBTYPE) %>% 
    summarise(count = n(), .groups = "drop") 
  
  p <- ggplot(mutation_summary) +
    geom_segment(aes(x=CHROM, xend=CHROM, y=0, yend=count), color="grey") +
    geom_point(aes(x=CHROM, y=count, color=SUBTYPE), size=3) +
    coord_flip()+
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",) +
    labs(title = "Subtype distribution on chromosomes")+
    facet_wrap(~SUBTYPE, ncol=1, scale="free_y")
  ggplotly(p)
}


# Distribúcia SNP typov
SNV_types <- function(vcf_tibble){
  vcf_tibble %<>% filter(TYPE == "SNP") %>%  
    count(SNP_TYPE = paste(REF, ALT, sep = ">"))
  
  p <- ggplot(vcf_tibble, aes(x = reorder(SNP_TYPE, -n), y = n, fill = SNP_TYPE)) +
    geom_bar(stat = "identity", color = "black") +
    theme_minimal() +
    labs(title = "SNP Mutations Distribution", y=NULL, x=NULL) +
    theme(legend.position = "none")+
    coord_flip()
  
  ggplotly(p)
}


# Heatmapa počtu mutácií na chromozómoch
mutation_heatmap <- function(vcf_tibble){
  vcf_tibble %<>% mutate(bucket = ceiling(POS/10000000)) %>% 
    group_by(CHROM, bucket) %>% 
    summarise(count = n(), .groups = "drop")
  
  p <- ggplot(vcf_tibble, aes(CHROM, bucket, fill= count)) + 
    geom_tile()+
    theme_minimal() +
    labs(title = "Mutation Heatmap") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))+
    scale_fill_gradient(low="skyblue", high="black")

  p
} 
