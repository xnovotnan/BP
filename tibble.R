library(dplyr)
library(tibble)
library(VariantAnnotation)
library(ggplot2)
library(BiocManager)
library(tidyverse)
library(magrittr)
library(purrr)

vcf_file_mini <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new_mini.vcf"
vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new.vcf"

#tibble bez grange
prepare_data <- function(vcf_path){
  vcf_lines <- readLines(vcf_path)
  vcf_lines <- vcf_lines[!grepl("^#", vcf_lines)]
  data_split <- strsplit(vcf_lines, "\t")
  
  vcf_tibble <- tibble(
    CHROM = sapply(data_split, function(x) x[1]),
    POS = as.integer(sapply(data_split, function(x) x[2])),
    REF = sapply(data_split, function(x) x[4]),
    ALT = sapply(data_split, function(x) x[5]),
    QUAL = as.numeric(sapply(data_split, function(x) x[6])),
    GT = sapply(data_split, function(x) strsplit(x[10], ":")[[1]][1]),
    AD = sapply(data_split, function(x) strsplit(x[10], ":")[[1]][2]),
    DP = as.integer(sapply(data_split, function(x) strsplit(x[10], ":")[[1]][3]))
  )  
  vcf_tibble %<>% filter(QUAL > 200) %>% filter(GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"))
  return(vcf_tibble)
}

a_200 <- prepare_data(vcf_file_mini)
dim(a_bez)
dim(a_30)
dim(a_100)
dim(a_200)

# tibble cez grange
prepare_data2 <- function(vcf_path){
  vcf <- readVcf(vcf_path, "t2t")
  vcf_tibble <- tibble(
    CHROM = as.character(seqnames(vcf)),
    POS = start(vcf),
    REF = as.character(ref(vcf)),
    ALT = rowRanges(vcf)$ALT %>% CharacterList(.) %>% unstrsplit(., sep = ","),
    QUAL = qual(vcf),
    GT = as.character(geno(vcf)$GT),
    AD = as.character(geno(vcf)$AD),
    DP = as.character(geno(vcf)$DP)
  )
  vcf_tibble %<>% filter(QUAL > 200) %>% filter(GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"))
  return(vcf_tibble)
}

b <- prepare_data2(vcf_file)
print(head(b))


#klasifikacia types
get_types <- function(vcf_tibble) {
  get_first_alt <- function(alt) { 
    alt_split <- strsplit(alt, ",")[[1]]  
    return(alt_split[1])
  }
  
  vcf_tibble %<>%
    mutate(TYPE = ifelse(
      nchar(REF) == 1 & nchar(sapply(ALT, get_first_alt)) == 1, 
      "SNP", "INDEL"
    ))
  return(vcf_tibble)
}

#klasifikacia subtypes
get_subtypes <- function(vcf_tibble) {
  get_first_alt <- function(alt) { 
    alt_split <- strsplit(alt, ",")[[1]]  
    return(alt_split[1])
  }
  purines <- c("A", "G")
  pyrimidines <- c("C", "T")

  vcf_tibble %<>%
    mutate(SUBTYPE = ifelse(
      nchar(REF) == 1 & nchar(sapply(ALT, get_first_alt)) == 1, 
      ifelse(
        (REF %in% purines & sapply(ALT, get_first_alt) %in% purines) | 
          (REF %in% pyrimidines & sapply(ALT, get_first_alt) %in% pyrimidines), 
        "Transition", "Transversion"), 
      ifelse(
        nchar(REF) > nchar(sapply(ALT, get_first_alt)),
        "Deletion", "Insertion"
      )))
  return(vcf_tibble)
}

a <- get_types(a)
a <- get_subtypes(a)
head(a)

b <- get_types(b)
b <- get_subtypes(b)
head(b)

#vypocet AF
get_frequency <- function(vcf_tibble) {
  ad <- sapply(strsplit(vcf_tibble$AD, ","), function(x){as.numeric(x[1])})
  vcf_tibble %<>%
    mutate(AF = map2_dbl(ad, vcf_tibble$DP, ~ {
    if (!is.null(.x) && !is.null(.y) && .y > 0) {
      1 - .x / .y
    } else {NA}}))
  
  return(vcf_tibble)
}


# distribucia mutacii na chromozomoch
get_distribution <- function(vcf_tibble, chromosome="all", subtypes=FALSE){
  if (chromosome != "all") vcf_tibble %<>% filter(CHROM == chromosome) 
  if (!subtypes) types <- vcf_tibble$TYPE
  else types <- vcf_tibble$SUBTYPE
  
  ggplot(vcf_tibble %>% 
           filter(chromosome == "all" | CHROM == chromosome), 
         aes(x = CHROM, fill = types)) +
    geom_bar(position = "stack", color = "black") +
    labs(title = "Variant Distribution on Chromosomes",
         x = "Chromosome",
         y = "Variant count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

get_distribution(b, chromosome="chr1")
get_distribution(b, subtypes=TRUE)
get_distribution(b, subtypes=FALSE)


#alelicka frekvencia 
plot_frequency <- function(vcf_tibble) {
  ggplot(vcf_tibble, aes(x=seq_along(AF), y = AF)) +
    geom_col(fill = "skyblue") +
    labs(title = "Allele Frequency",
         x = "Mutation",
         y = "Frequency") +
    theme_minimal()
  
}

a <- get_frequency(a)
plot_frequency(a)


# distribucia mutacii celkovo
mut_summary <- function(vcf_tibble, subtypes=FALSE){
  if(subtypes) vcf_tibble %<>% mutate(TYPE = vcf_tibble$SUBTYPE)
  ggplot(vcf_tibble, aes(x=TYPE, fill = TYPE)) +
    geom_bar(color = "black") + 
    labs(title = "Mutation Types Distribution",
         x = "Mutation Type",
         y = "Frequency") +
    theme_minimal()+
    coord_flip()
}

mut_summary(a)
mut_summary(a, subtypes=TRUE)


# SNV typy
SNV_types <- function(vcf_tibble){
  vcf_tibble %<>% filter(TYPE == "SNP") %>%  
    mutate(SNP_TYPE = map2_chr(REF, ALT, ~{
      paste(.x, ">", .y)
    }))
  
  ggplot(vcf_tibble, aes(x = SNP_TYPE, fill = SNP_TYPE)) +
    geom_bar(color = "black") + 
    labs(title = "Distribution of SNP Types",
         x = "SNP Type",
         y = "Frequency") +
    theme_minimal() +
    coord_flip()
}

SNV_types(a)
head(a)


# Boxplotymutacii na chromozomoch
mut_summary_boxplot <- function(vcf_tibble, subtypes = FALSE) {
  if (subtypes) vcf_tibble <- vcf_tibble %>% mutate(TYPE = SUBTYPE)
  
  mutation_counts <- vcf_tibble %>%
    group_by(CHROM, TYPE) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(CHROM) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup() %>%
    group_by(TYPE)
  
  print(mutation_counts)
  
  ggplot(mutation_counts, aes(x = TYPE, y = percent, fill = TYPE)) +
    geom_boxplot() +
    labs(title = "Boxplot of Mutation Percentages by Type",
         x = "Mutation Type",
         y = "Percentage of Mutations (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_flip()
}

mut_summary_boxplot(b)

