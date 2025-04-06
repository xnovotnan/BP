library(tidyverse)
library(magrittr)
library(ggridges)
library(scales)
library(utils)
library(vroom)
options(scipen = 999)

# Preprocess vcf data
process_vcf <- function(vcf_file){
  name <- basename(vcf_file)
  family <- stringr::str_extract(name, "(?<=\\.)\\d+")
  vcf_tibble <- vroom::vroom(
    file = vcf_file,
    delim = "\t",
    comment = "##"
  )

  vcf_tibble %<>%
    rename("VALUES" = last(colnames(.)), "CHROM" = first(colnames(.))) %>%
    mutate(
      GT = str_split(VALUES, ":") %>% map_chr(~ .x[1]),
      AD = str_split(VALUES, ":") %>% map_chr(~ .x[2]) %>% str_split(",") %>% map_int(~ as.integer(.x[1])),
      DP = str_split(VALUES, ":") %>% map_int(~ as.integer(.x[3])),
      TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"SNP", "INDEL"),
      AF = round(AD/DP, 2)) %>%
    filter(QUAL > 200,
           GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"),
           DP>0,
           CHROM != "chrM")

  purines <- c("A", "G")
  pyrimidines <- c("C", "T")
  vcf_tibble %<>% mutate(
    SUBTYPE = ifelse(TYPE == "SNP",
                     ifelse((REF %in% purines & ALT %in% purines) | (REF %in% pyrimidines & ALT %in% pyrimidines),
                            "Transition",
                            "Transversion"),
                     ifelse(nchar(REF) > nchar(ALT),
                            "Deletion",
                            "Insertion"))
    )

  vcf_tibble$CHROM <- factor(vcf_tibble$CHROM, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
  vcf_tibble %<>% select(c("CHROM", "POS", "REF", "ALT", "TYPE", "SUBTYPE", "QUAL", "GT", "AD", "DP", "AF"))
  tibble(
    name = name,
    family = family,
    data = list(vcf_tibble)
  )
}
prepare_vcf_files <- function(folder_path){
  files <- list.files(folder_path, pattern = "\\.vcf.gz$", full.names = TRUE)
  vcf_data <- purrr::map_dfr(files, process_vcf)
  vcf_data
}
# b <- prepare_vcf_files("/Users/macbook/Documents/BP/data/testik.vcf")

# Mutation Counts
num_of_mutation <- function(vcf_data){
  vcf_data %<>% mutate(mutation_count = map_int(data, nrow))

  p <- ggplot(vcf_data, aes(x = mutation_count, y = reorder(name, mutation_count), fill = name)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title =  element_blank(),
         x = element_blank(),
         y= element_blank()) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12)) +
    scale_fill_brewer(palette = "Set3")
  p
}
mutation_heatmap <- function(vcf_data){
  vcf_data %<>%
    mutate(data = map2(data, name, ~mutate(.x, SAMPLE = .y))) %>%
    unnest(data) %>%
    group_by(SAMPLE, CHROM) %>%
    summarise(count = n(),
              max_pos = max(POS),
              .groups = "drop") %>%
    mutate(percentage = count / max_pos) %>%
    select(SAMPLE, CHROM, percentage)
  
  p <- ggplot(vcf_data, aes(CHROM, SAMPLE, fill= percentage)) + 
    geom_tile() +
    labs(title = "Mutation Counts across Chromosomes",
         x = element_blank(),
         y = element_blank()) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12, hjust = 1, angle = 90),
          axis.text.y = element_text(size = 12),
          legend.key.width = unit(1, "cm"))
  p
}
mutation_types_distribution <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    mutate(data = map2(data, name, ~mutate(.x, SAMPLE = .y))) %>%
    unnest(data) %>%
    select(SAMPLE, TYPE) %>%
    group_by(SAMPLE, TYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(SAMPLE) %>%
    mutate(percentage = count * 100/ sum(count))
  mutationCounts$TYPE <- with(mutationCounts, reorder(TYPE, percentage, FUN = sum))
  
  p <- ggplot(mutationCounts, aes(x = SAMPLE, y = percentage, fill = TYPE)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, title = "Variant Distribution on Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank())
  p
}
mutation_subtypes_distribution <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    mutate(data = map2(data, name, ~mutate(.x, SAMPLE = .y))) %>%
    unnest(data) %>%
    select(SAMPLE, SUBTYPE) %>%
    group_by(SAMPLE, SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(SAMPLE) %>%
    mutate(percentage = count * 100/ sum(count))
  mutationCounts$SUBTYPE <- with(mutationCounts, reorder(SUBTYPE, percentage, FUN = sum))
  
  p <- ggplot(mutationCounts, aes(x = SAMPLE, y = percentage, fill = SUBTYPE)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, title = "Variant Distribution on Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank())
  p
}

# SNP Analysis
snp_class_comparison <- function(vcf_data){
  vcf_data %<>% mutate(data = map2(data, name, ~mutate(.x, SAMPLE = .y))) %>%
    unnest(data) %>%
    filter(TYPE == "SNP") %>%
    select(SAMPLE, REF, ALT) %>%
    group_by(SAMPLE) %>%
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  vcf_data$SNP_TYPE <- with(vcf_data, reorder(SNP_TYPE, percentage, FUN = sum))
  
  p <- ggplot(vcf_data, aes(x = percentage, y = SAMPLE, fill = SNP_TYPE)) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage", 
         y = element_blank(), 
         title= "SNP Types Distribution on Samples") +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "right")
  p
}
transversion_transitions <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    mutate(data = map2(data, name, ~mutate(.x, SAMPLE = .y))) %>%
    unnest(data) %>%
    filter(TYPE == "SNP") %>%
    select(SAMPLE, SUBTYPE) %>%
    group_by(SAMPLE, SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(SAMPLE) %>%
    mutate(percentage = count * 100/ sum(count))
  mutationCounts$SUBTYPE <- with(mutationCounts, reorder(SUBTYPE, percentage, FUN = sum))
  
  p <- ggplot(mutationCounts, aes(x = SAMPLE, y = percentage, fill = SUBTYPE)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, title = "Variant Distribution on Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank())
  p
}

# INDEL Analysis
insertion_deletions <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    mutate(data = map2(data, name, ~mutate(.x, SAMPLE = .y))) %>%
    unnest(data) %>%
    filter(TYPE == "INDEL") %>%
    select(SAMPLE, SUBTYPE) %>%
    group_by(SAMPLE, SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(SAMPLE) %>%
    mutate(percentage = count * 100/ sum(count))
  mutationCounts$SUBTYPE <- with(mutationCounts, reorder(SUBTYPE, percentage, FUN = sum))
  
  p <- ggplot(mutationCounts, aes(x = SAMPLE, y = percentage, fill = SUBTYPE)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = NULL, y = NULL, title = "Variant Distribution on Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_minimal()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_blank())
  p
}
indel_len_boxplot <- function(vcf_data){
  vcf_data %<>% 
    mutate(data = map2(data, name, ~mutate(.x, sample = .y))) %>%
    unnest(data) %>%
    filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
    
  p <- ggplot(vcf_data, aes(x = length, y = sample, fill = sample)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = "Boxplot of INDEL Length per Sample",  x = "INDEL Length", y = element_blank()) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}

# Quality Analysis
quality_boxplot_per_fam <- function(vcf_data){
  plot_data <- vcf_data %>%
    mutate(data = map2(data, name, ~mutate(.x, sample = .y))) %>%
    unnest(data) %>%
    select(sample, family, QUAL)
  
  p <- ggplot(plot_data, aes(x = QUAL, y = family, fill=family)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = "Boxplot of Quality per Family",  x = "Quality", y = element_blank()) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}
quality_boxplot <- function(vcf_data){
  vcf_data %<>%
    mutate(data = map2(data, name, ~mutate(.x, sample = .y))) %>%
    unnest(data) %>%
    select(sample, QUAL)
  
  p <- ggplot(vcf_data, aes(x = QUAL, y = sample, fill = sample)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = "Boxplot of QUAL per Sample",  x = "Quality", y = element_blank()) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}


# Allele Frequency Analysis
frequency_ridges <- function(vcf_data){
  vcf_data %<>%
    mutate(data = map2(data, name, ~mutate(.x, sample = .y))) %>%
    unnest(data) %>%
    select(sample, AF)
  
  p <- ggplot(vcf_data, aes(x=AF, y=sample, fill=sample)) +
    geom_density_ridges() +
    labs(title = "Allele Frequency across the chromosomes",
         x= element_blank(),
         y=element_blank()) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  p
}

# Read Depth Analysis
read_depth <- function(vcf_data){
  vcf_data %<>%
    mutate(data = map2(data, name, ~mutate(.x, sample = .y))) %>%
    unnest(data) %>%
    select(sample, DP)
  
  p <- ggplot(vcf_data, aes(x=DP, y=sample, fill= sample)) +
    geom_violin(color="black") +
    labs(title = "Read depth across chromosomes",
         x = element_blank(),
         y = element_blank()) +
    theme_minimal() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size=12),
          axis.text.y = element_text(size = 12))
  p
}




