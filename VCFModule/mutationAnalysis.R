library(tidyverse)
library(magrittr)
library(ggridges)
library(scales)
library(ggrepel)
library(plotly)
options(scipen = 999)

# The prepare_data function reads the VCF file using the vroom package. The 
# relevant data is extracted from the file. The mutations are classified into 
# mutation types (SNP or INDEL) and subtypes (transitions or transversions, 
# insertion, deletion, or structural variant) based on the REF and ALT field.
# The data is then filtered based on QUAL, GT, and DP thresholds. The AF is 
# calculated from the AD and DP values.
prepare_data <- function(vcf_file){
  incProgress(0.1, detail = "Reading file...")
  vcf_tibble <- vroom::vroom(
    file = vcf_file,
    delim = "\t",
    comment = "##"
  )
  incProgress(0.3, detail = "Processing data...")
  vcf_tibble %<>% 
    rename("VALUES" = last(colnames(.)), "CHROM" = first(colnames(.))) %>%
    mutate(
      GT = str_split(VALUES, ":") %>% map_chr(~ .x[1]),
      AD = str_split(VALUES, ":") %>% map_chr(~ .x[2]) %>% str_split(",") %>% map_int(~ as.integer(.x[1])),
      DP = str_split(VALUES, ":") %>% map_int(~ as.integer(.x[3])),
      TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"SNP", "INDEL")) %>%
    filter(QUAL > 50, 
           GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"),
           DP>0,
           CHROM != "chrM")
  incProgress(0.7, detail = "Categorizing mutations...")
  purines <- c("A", "G")
  pyrimidines <- c("C", "T")
  vcf_tibble %<>% 
    mutate(
      SUBTYPE = ifelse(TYPE == "SNP", 
                       ifelse((REF %in% purines & ALT %in% purines) | 
                                (REF %in% pyrimidines & ALT %in% pyrimidines), "Transition", "Transversion"), 
                       ifelse(abs(nchar(REF) - nchar(ALT)) > 50, "Structural Variant",
                         ifelse(nchar(REF) > nchar(ALT),"Deletion", "Insertion"))
                       ),
      AF = round(AD/DP, 2)
    )
  incProgress(0.9, detail = "Completing preprocessing...")
  vcf_tibble$CHROM <- factor(vcf_tibble$CHROM, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
  vcf_tibble %<>% select(c("CHROM", "POS", "REF", "ALT", "TYPE", "SUBTYPE", "QUAL", "GT", "AD", "DP", "AF"))
  vcf_tibble
}

# Mutation Counts Section
mutation_count <- function(vcf_tibble){
  sprintf("Mutation Count: %s", label_comma()(dim(vcf_tibble)[1]))
}
mutation_donut <- function(vcf_tibble, subtypes=FALSE){
  if (subtypes) vcf_tibble$TYPE <- vcf_tibble$SUBTYPE
  # calculates percentual distribution of mutation types
  mutationCounts <- vcf_tibble %>%
    group_by(TYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(count * 100 / sum(count), 2))
  
  # calculates labels positions
  mutationCounts %<>%
    mutate(csum = rev(cumsum(rev(percentage))), 
           pos = percentage / 2 + lead(csum, 1),
           pos = if_else(is.na(pos), percentage / 2, pos))
  
  p <- ggplot(mutationCounts, aes(x = "" , y = percentage, fill = TYPE)) +
    geom_col(width = 1, color = 1, show.legend = FALSE) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Set3") +
    geom_label_repel(data = mutationCounts,
                     aes(y = pos, label = paste0(TYPE, "s\n", percentage, "%\n", label_comma()(count))),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    labs(title = "Mutation Type Distribution")
  p
}
mutation_distribution <- function(vcf_tibble, subtypes=FALSE){
  if (subtypes) vcf_tibble$TYPE <- vcf_tibble$SUBTYPE
  # calculates percentual distribution of mutation types per chromosome
  vcf_tibble %<>% 
    group_by(CHROM, TYPE) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    group_by(CHROM) %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  vcf_tibble$TYPE <- with(vcf_tibble, reorder(TYPE, percentage))
  
  p <- ggplot(vcf_tibble, aes(x = CHROM, 
                              y = percentage, 
                              fill = TYPE,
                              text = paste0("Type: ", TYPE, "\n",
                                            "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = element_blank(), y = "Percentage", title = "Mutation Type Distribution Across Chromosomes (%)") +
    coord_flip() +
    scale_fill_brewer(palette = "Set3")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
  p
}

# SNP Analysis Section
snp_value <- function(vcf_tibble){
  # calculates SNP count
  values <- vcf_tibble %>% 
    group_by(TYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2)) %>%
    filter(TYPE == "SNP")
  sprintf("SNP Count: %s (%s %%)", label_comma()(values$count), values$percentage)
}
snp_types_donut <- function(vcf_tibble){
  # calculates percentual distribution of SNP subtypes
  vcf_tibble %<>% filter(TYPE == "SNP")
  mutationCounts <- vcf_tibble %>%
    group_by(SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  
  # calculates labels positions
  mutationCounts <- mutationCounts %>%
    mutate(csum = rev(cumsum(rev(percentage))), 
           pos = percentage / 2 + lead(csum, 1),
           pos = if_else(is.na(pos), percentage / 2, pos))
  
  p <- ggplot(mutationCounts, aes(x = "" , y = percentage, fill = SUBTYPE)) +
    geom_col(width = 1, color = 1, show.legend = FALSE) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Set3") +
    geom_label_repel(data = mutationCounts,
                     aes(y = pos, label = paste0(SUBTYPE, "s\n", percentage, "%\n", label_comma()(count))),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    labs(title = "SNP Type Distribution")
  p
}
snp_class_stacked <- function(vcf_tibble){
  # calculates percentual distribution of SNP substitution types (for example A>C, C>G, ...)
  vcf_tibble %<>% 
    filter(TYPE == "SNP") %>% 
    group_by(CHROM) %>%
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2))
  vcf_tibble$SNP_TYPE <- with(vcf_tibble, reorder(SNP_TYPE, percentage))
  
  p <- ggplot(vcf_tibble, aes(x = percentage, y = CHROM, fill = SNP_TYPE,
                              text = paste0("Type: ", SNP_TYPE, "\n",
                                            "Chromosome: ", CHROM, "\n",
                                            "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage",
         y = element_blank(), 
         title= "SNP Substitution Type Distribution Across Chromosomes (%)") +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
  p
}
snp_class_boxplot <- function(vcf_tibble){
  # calculates percentual distribution of SNP substitution types per chromosome 
  # and classifies them as transversions or trasitions
  vcf_tibble %<>% 
    filter(TYPE == "SNP") %>% 
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>%
    mutate(percentage = round(n *100/ sum(n), 2),
           subtype = ifelse(SNP_TYPE %in% c("A>G", "G>A", "C>T", "T>C"), 
                            "Transition", "Transversion"))
  vcf_tibble$SNP_TYPE <- with(vcf_tibble, reorder(SNP_TYPE, percentage))
  
  p <- ggplot(vcf_tibble, 
              aes(x = reorder(SNP_TYPE, -percentage), 
                  y = percentage, 
                  color = subtype, 
                  fill = SNP_TYPE)) +
    geom_boxplot() +
    labs(title = "Distribution of SNP Substitution Types (%)",
         x = "Base Substitution",
         y = "Percentage") +
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          legend.title = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x =  element_text(size = 14),
          axis.title.y =  element_text(size = 14)) +
    scale_fill_brewer(palette = "Set3")
  p
}
snp_class_barplot <- function(vcf_tibble){
  # calculates percentual distribution of SNP substitution types
  vcf_tibble %<>% 
    filter(TYPE == "SNP") %>% 
    count(SNP_TYPE = paste(REF, ALT, sep = ">")) %>% 
    mutate(percentage = round(n *100/ sum(n), 2),
           label = paste0(SNP_TYPE, "\n", percentage, "%"))
  vcf_tibble$SNP_TYPE <- with(vcf_tibble, reorder(SNP_TYPE, percentage))
  
  p <- ggplot(vcf_tibble, aes(x = percentage, y = reorder(SNP_TYPE, n), fill = SNP_TYPE)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    labs(title = "SNP Substitution Type Distribution (%)") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    geom_text(aes(x = percentage/2, label = label), size = 4, color = "black")+
    scale_fill_brewer(palette = "Set3")
  p
}

# INDEL Analysis Section
indel_values <- function(vcf_tibble){
  values <- vcf_tibble %>% 
    group_by(TYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2)) %>%
    filter(TYPE == "INDEL")
  sprintf("INDEL Count: %s (%s %%)", label_comma()(values$count), values$percentage)
}
indel_types <- function(vcf_tibble){
  # calculates percentual distribution of INDEL subtypes
  vcf_tibble %<>% 
    filter(TYPE == "INDEL")
  mutationCounts <- vcf_tibble %>%
    group_by(SUBTYPE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(percentage = round(count * 100/ sum(count),2))
  
  # calculates labels positions
  mutationCounts <- mutationCounts %>%
    mutate(csum = rev(cumsum(rev(percentage))), 
           pos = percentage / 2 + lead(csum, 1),
           pos = if_else(is.na(pos), percentage / 2, pos))
  
  p <- ggplot(mutationCounts, aes(x = "" , y = percentage, fill = SUBTYPE)) +
    geom_col(width = 1, color = 1, show.legend = FALSE) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Set3") +
    geom_label_repel(data = mutationCounts,
                     aes(y = pos, label = paste0(SUBTYPE, "s\n", percentage, "%\n", label_comma()(count))),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) +
    labs(title = "Insertion-Deletion Distribution")
  p
}
indel_stacked <- function(vcf_tibble){
  # calculates percentual distribution of INDEL subtypes per chromosome
  vcf_tibble %<>% 
    filter(TYPE == "INDEL") %>% 
    group_by(CHROM, SUBTYPE) %>%
    summarise(count = n()) %>%
    mutate(percentage = round(count *100/ sum(count), 2))
  vcf_tibble$SUBTYPE <- with(vcf_tibble, reorder(SUBTYPE, -percentage))
  
  p <- ggplot(vcf_tibble, aes(x = percentage, y = CHROM, fill = SUBTYPE,
                              text = paste0("Type: ", SUBTYPE, "\n",
                                            "Chromosome: ", CHROM, "\n",
                                            "Percentage: ", percentage, "%"))) +
    geom_bar(stat = "identity", position = "stack") + 
    labs(x = "Percentage", 
         y = element_blank(), 
         title= "Insertions and Deletions Across Chromosomes (%)") +
    theme_classic() +
    scale_fill_brewer(palette = "Set3") +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_blank())
  p
}
indel_length_avg <- function(vcf_tibble){
  # calculates mean INDEL size
  vcf_tibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  sprintf("Mean INDEL size: %s", round(mean(vcf_tibble$length),2))
}
indel_length_med <- function(vcf_tibble){
  # calculates median INDEL size
  vcf_tibble %<>% filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF)))
  sprintf("Median size length: %s", median(vcf_tibble$length))
}
indel_length <- function(vcf_tibble){
  # filters out strucutral variants and calculates INDEL size
  vcf_tibble %<>%
    filter(TYPE == "INDEL", SUBTYPE != "Structural Variant") %>%
    mutate(length = nchar(ALT) - nchar(REF))

  p <- ggplot(vcf_tibble, aes(x = length, fill=SUBTYPE)) +
    geom_histogram(binwidth = 3,
                   position = "dodge",
                   alpha = 0.6) +
    labs(x = "Length (bp)",
         y = "Mutation Count",
         title = "Distribution of INDEL Sizes") +
    theme_classic() +
    scale_y_continuous(labels = comma) +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
  p
}
indel_length_boxplot <- function(vcf_tibble){
  # calculates INDEL sizes per subtype
  vcf_tibble %<>% 
    filter(TYPE == "INDEL") %>%
    mutate(length = abs(nchar(ALT) - nchar(REF))) %>%
    group_by(CHROM, SUBTYPE)
  
  p <- ggplot(vcf_tibble, aes(x = SUBTYPE, y = length, fill = SUBTYPE)) +
    geom_boxplot(show.legend = FALSE) +
    labs(title= "INDEL Length Distribution by Type",
         x = "Type", 
         y = "Length (bp)") +
    theme_classic() +
    scale_fill_brewer(palette = "Set3")+
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.title.x =  element_text(size = 14),
          axis.title.y =  element_text(size = 14))
  p
}

# Quality Section
quality_avg <- function(vcf_tibble){
  sprintf("Mean Quality: %s", round(mean(vcf_tibble$QUAL),2))
}
quality_med <- function(vcf_tibble){
  sprintf("Median Quality: %s", median(vcf_tibble$QUAL))
}
quality_bar <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=seq_along(QUAL), y = QUAL)) +
    geom_hex() +
    labs(title = "Quality Across the Genome", 
         x = "Position", 
         y="Quality") +
    theme_classic() +
    scale_x_continuous(labels = comma) +
    scale_fill_gradientn(colors = c("lightgreen", "green4")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_blank())
  p
}
quality_on_chroms <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=CHROM, y=QUAL)) +
    geom_violin(color="orange") +
    labs(title = "Quality Across Chromosomes",
         x = element_blank(),
         y = "Quality") +
    theme_classic() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size=14),
          axis.text.y = element_text(size = 14),
          axis.title.y =  element_text(size = 14))
  p
}

# Allele Frequency Section
allele_freq_avg <- function(vcf_tibble){
  sprintf("Mean Allele Frequency: %s", round(mean(vcf_tibble$AF),2))
}
allele_freq_med <- function(vcf_tibble){
  sprintf("Median Allele Frequency: %s", median(vcf_tibble$AF))
}
allele_freq_hexbin <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=seq_along(AF), y = AF)) +
    geom_hex() +
    labs(title = "Allele Frequency Across the Genome", 
         x = "Position", 
         y="Allele Frequency") +
    theme_classic() +
    scale_x_continuous(labels = comma) +
    scale_fill_gradientn(colors = c("skyblue1", "blue4")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_blank())
  p
}
allele_freq_on_chroms <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=AF, y=CHROM)) +
    geom_density_ridges(fill="lightpink") +
    labs(title = "Allele Frequency Across Chromosomes",
         x= "Allele Frequency",
         y=element_blank()) +
    theme_classic() +
    theme(legend.position = "none", 
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x =  element_text(size = 14))
  p
}

# Read Depth Section
read_depth_avg <- function(vcf_tibble){
  sprintf("Mean Read Depth: %s", round(mean(vcf_tibble$DP),2))
}
read_depth_med <- function(vcf_tibble){
  sprintf("Median Read Depth: %s", median(vcf_tibble$DP))
}
read_depth_density <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=seq_along(DP), y = DP)) +
    geom_hex() +
    labs(title = "Read Depth Across the Genome", 
         x = "Position", 
         y="Read Depth") +
    theme_classic() +
    scale_x_continuous(labels = comma) +
    scale_fill_gradientn(colors = c("lightgreen", "green4")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10),
          legend.title = element_blank())
  p
}
read_depth_on_chroms <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=CHROM, y=DP)) +
    geom_violin(color="orange") +
    labs(title = "Read depth Across Chromosomes",
         x = element_blank(),
         y = "Read Depth") +
    theme_classic() +
    theme(legend.position = "none",  
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size=14),
          axis.text.y = element_text(size = 14),
          axis.title.y =  element_text(size = 14))
  p
}
