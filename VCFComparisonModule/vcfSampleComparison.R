library(tidyverse)
library(magrittr)
library(ggridges)
library(scales)
library(utils)
library(vroom)
options(scipen = 999)
library(eulerr)
library(grid)
library(gridExtra)

# Preprocess vcf data
process_vcf <- function(vcf_file){
  incProgress(0.1, detail = "Processing data...")
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
    name = sub("\\.vcf\\.gz$", "", name),
    family = family,
    data = list(vcf_tibble)
  )
}
prepare_vcf_files <- function(folder_path){
  files <- list.files(folder_path, pattern = "\\.vcf.gz$", full.names = TRUE)
  vcf_data <- purrr::map_dfr(files, process_vcf)
  incProgress(0.9, detail = "Completing preprocessing...")
  vcf_data
}

# Mutation Counts
num_of_mutation <- function(vcf_data){
  vcf_data %<>% mutate(mutation_count = map_int(data, nrow),
                       label = paste0(name, " - \n", label_comma()(mutation_count)))

  p <- ggplot(vcf_data, aes(x = mutation_count, y = reorder(name, mutation_count), fill = name)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "Mutation Counts per Sample") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))+
    geom_text(aes(x = mutation_count/2, label = label), size = 5, color = "black")+
    scale_fill_brewer(palette = "Set3")
  p
}
mutation_heatmap <- function(vcf_data){
  vcf_data %<>%
    unnest(data) %>%
    group_by(name, CHROM) %>%
    summarise(count = n(), .groups = "drop") %>%
    select(name, CHROM, count)
  
  p <- ggplot(vcf_data, aes(CHROM, name, fill= count,
                            text = paste0("Sample: ", name, "\n",
                                          "Chromosome: ", CHROM, "\n",
                                          "Count: ", count))) + 
    geom_tile() +
    labs(title = "Mutation Counts across Chromosomes",
         x = element_blank(),
         y = element_blank()) +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10, hjust = 1, angle = 45),
          axis.text.y = element_text(size = 10),
          legend.key.width = unit(1, "cm"))+
    scale_fill_gradient(low = "lightgreen", high = "orange")
  p
}
mutation_types_distribution <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    unnest(data) %>%
    select(name, TYPE) %>%
    group_by(name, TYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(name) %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  mutationCounts$TYPE <- with(mutationCounts, reorder(TYPE, percentage))
  
  p <- ggplot(mutationCounts, aes(x = name, y = percentage, fill = TYPE,
                                  text = paste0("Sample: ", name, "\n",
                                                "Type: ", TYPE, "\n",
                                                "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = element_blank(), y = "Percentage", title = "Distribution of Mutation Types Across Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank())
  p
}
mutation_subtypes_distribution <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    unnest(data) %>%
    select(name, SUBTYPE) %>%
    group_by(name, SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(name) %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  mutationCounts$SUBTYPE <- with(mutationCounts, reorder(SUBTYPE, percentage))
  
  p <- ggplot(mutationCounts, aes(x = name, y = percentage, fill = SUBTYPE,
                                  text = paste0("Sample: ", name, "\n",
                                                "Type: ", SUBTYPE, "\n",
                                                "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = element_blank(), y = "Percentage", title = "Distribution of Mutation Subtypes Across Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank())
  p
}

# SNP Analysis
snp_class_comparison <- function(vcf_data){
  vcf_data %<>% 
    unnest(data) %>%
    filter(TYPE == "SNP") %>%
    select(name, REF, ALT) %>%
    group_by(name) %>%
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  vcf_data$SNP_TYPE <- with(vcf_data, reorder(SNP_TYPE, percentage))

    p <- ggplot(vcf_data, aes(x = percentage, 
                              y = name, 
                              fill = SNP_TYPE, 
                              text = paste0("Sample: ", name, "<br>",
                                            "Type: ", SNP_TYPE, "<br>",
                                            "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage", 
         y = element_blank(), 
         title= "SNP Substitution Type Distribution on Samples") +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank(),
          legend.position = "right")
  p
}
transversion_transitions <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    unnest(data) %>%
    filter(TYPE == "SNP") %>%
    select(name, SUBTYPE) %>%
    group_by(name, SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(name) %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  mutationCounts$SUBTYPE <- with(mutationCounts, reorder(SUBTYPE, percentage, FUN = sum))
  
  p <- ggplot(mutationCounts, aes(x = name, y = percentage, fill = SUBTYPE,
                                  text = paste0("Sample: ", name, "\n",
                                                "Type: ", SUBTYPE, "\n",
                                                "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = element_blank(), y = "Percentage", title = "Transversions-Transitions Distribution on Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank())
  p
}

# INDEL Analysis
insertion_deletions <- function(vcf_data){
  mutationCounts <- vcf_data %>%
    unnest(data) %>%
    filter(TYPE == "INDEL") %>%
    select(name, SUBTYPE) %>%
    group_by(name, SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(name) %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  mutationCounts$SUBTYPE <- with(mutationCounts, reorder(SUBTYPE, percentage, FUN = sum))
  
  p <- ggplot(mutationCounts, aes(x = name, y = percentage, fill = SUBTYPE,
                                  text = paste0("Sample: ", name, "\n",
                                                "Type: ", SUBTYPE, "\n",
                                                "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = element_blank(), y = "Percentage", title = "Insertion-Deletion Distribution on Samples (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank())
  p
}
indel_len_boxplot <- function(vcf_data){
  vcf_data %<>% 
    unnest(data) %>%
    filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
    
  p <- ggplot(vcf_data, aes(x = length, y = name, fill = name)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = "INDEL Size per Sample",  x = "Length (bp)", y = element_blank()) +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x =  element_text(size = 14))
  p
}

# Quality Analysis
quality_boxplot <- function(vcf_data){
  vcf_data %<>%
    unnest(data) %>%
    select(name, QUAL)
  
  p <- ggplot(vcf_data, aes(x = QUAL, y = name, fill = name)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title = "Quality per Sample",  x = "Quality", y = element_blank()) +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x =  element_text(size = 14))+
    scale_x_continuous(labels = comma)
  
  p
}

# Allele Frequency Analysis
frequency_ridges <- function(vcf_data){
  vcf_data %<>%
    unnest(data) %>%
    select(name, AF)
  
  p <- ggplot(vcf_data, aes(x=AF, y=name, fill=name)) +
    geom_density_ridges() +
    labs(title = "Allele Frequency Across Samples",
         x= "Allele Frequency",
         y=element_blank()) +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x =  element_text(size = 14))
  p
}

# Read Depth Analysis
read_depth <- function(vcf_data){
  vcf_data %<>%
    unnest(data) %>%
    select(name, DP)
  
  p <- ggplot(vcf_data, aes(x=DP, y=name, fill= name)) +
    geom_violin(color="black") +
    labs(title = "Read depth Across Samples",
         x = "Read depth",
         y = element_blank()) +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size = 14),
          axis.title.x =  element_text(size = 14))
  p
}

# VENN DIAGRAM 2 FILES
venn_diagram <- function(file1, file2, data){
  incProgress(0.2, detail = "Processing first file...")
  data1 <- data %>% filter(name == file1) %>%
    pull(data) %>%
    .[[1]] %>%
    mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))%>%
    select(ID)
  
  incProgress(0.2, detail = "Processing second file...")
  data2 <- data %>% filter(name == file2) %>%
    pull(data) %>%
    .[[1]]%>%
    mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))%>%
    select(ID)
  
  incProgress(0.15, detail = "Computing overlaps...")
  n_only1 <-  setdiff(data1, data2) %>% nrow()
  n_only2 <- setdiff(data2, data1) %>% nrow()
  n_common <- intersect(data1, data2) %>% nrow()

  incProgress(0.2, detail = "Comuting overlaps...")
  fit <- euler(c(
    "VCF1" = n_only1,
    "VCF2" = n_only2,
    "VCF1&VCF2" = n_common
  ))
  incProgress(0.2, detail = "Generating venn diagram...")
  
  venn_plot <- plot(fit,
                    fills = c("skyblue", "orange"),
                    labels = TRUE,
                    quantities = TRUE
  )
  venn_plot
}

# VENN DIAGRAM 3 FILES
venn_diagram_3files <- function(file1, file2, file3, data) {
  incProgress(0.15, detail = "Processing first file...")
  data1 <- data %>% filter(name == file1) %>%
    pull(data) %>%
    .[[1]] %>%
    mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))%>%
    select(ID)
  
  incProgress(0.15, detail = "Processing second file...")
  data2 <- data %>% filter(name == file2) %>%
    pull(data) %>%
    .[[1]]%>%
    mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))%>%
    select(ID)

  incProgress(0.15, detail = "Processing third file...")
  data3 <- data %>% filter(name == file3) %>%
    pull(data) %>%
    .[[1]]%>%
    mutate(ID = paste0(CHROM, "_", POS, "_", REF, "_", ALT))%>%
    select(ID)
  
  n_only1 <-  setdiff(data1, union(data2, data3)) %>% nrow()
  n_only2 <- setdiff(data2, union(data1, data3)) %>% nrow()
  n_only3 <- setdiff(data3, union(data1, data2)) %>% nrow()
  n_common12 <- intersect(data1, data2) %>% setdiff(data3) %>% nrow()
  n_common13 <- intersect(data1, data3) %>% setdiff(data2) %>% nrow()
  n_common23 <- intersect(data2, data3) %>% setdiff(data1) %>% nrow()
  n_common123 <- Reduce(intersect, list(data1, data2, data3)) %>% nrow()
  
  incProgress(0.15, detail = "Computing overlaps...")
  fit <- euler(c(
    "VCF1" = n_only1,
    "VCF2" = n_only2,
    "VCF3" = n_only3,
    "VCF1&VCF2" = n_common12,
    "VCF1&VCF3" = n_common13,
    "VCF2&VCF3" = n_common23,
    "VCF1&VCF2&VCF3" = n_common123
  ))
  incProgress(0.3, detail = "Generating venn diagram...")
  venn_plot <- plot(fit,
                    fills = c("skyblue", "orange", "green3"),
                    labels = TRUE,
                    quantities = TRUE
  )
  
  venn_plot
}



