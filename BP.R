library(ggplot2)
library(BiocManager)
library(VariantAnnotation)
library(tidyverse)
library(magrittr)
library(dplyr)
library(purrr)

vcf_file_mini <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new_mini.vcf"
vcf_mini <- readVcf(vcf_file_mini, "t2t")
mcols(vcf_mini)$alt <- rowRanges(vcf_mini)$ALT %>% CharacterList(.) %>% unstrsplit(., sep = ",")
#alt(vcf_mini)[[22]] 

vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new.vcf"
vcf <- readVcf(vcf_file, "t2t")
mcols(vcf)$alt <- rowRanges(vcf)$ALT %>% CharacterList(.) %>% unstrsplit(., sep = ",")
rowRanges(vcf_mini)


#klasifikacia types
get_types <- function(vcf_data){
  get_first_alt <- function(alt) { alt[[1]] }
   mut_types <- ifelse(nchar(as.character(ref(vcf_data))) == 1 & 
           nchar(sapply(mcols(vcf_data)$alt, get_first_alt)) == 1, 
         "SNP","INDEL")
   return(as.character(mut_types))
}

#klasifikacia subtypes
get_subtypes <- function(vcf_data){
  get_first_alt <- function(alt) { alt[[1]] }
  purines <- c("A", "G")
  pyrimidines <- c("C", "T")
  mut_subtypes <- ifelse(nchar(as.character(ref(vcf_data))) == 1 & nchar(sapply(mcols(vcf_data)$alt, get_first_alt)) == 1, 
                        ifelse(((ref(vcf_data) %in% purines & sapply(mcols(vcf_data)$alt, get_first_alt) %in% purines) | 
                                  (ref(vcf_data) %in% pyrimidines & sapply(mcols(vcf_data)$alt, get_first_alt) %in% pyrimidines)), 
                               "Transition", "Transversion"), 
                        ifelse(nchar(as.character(ref(vcf_data))) > nchar(sapply(mcols(vcf_data)$alt, get_first_alt)),
                               "Deletion", "Insertion"))
  return(as.character(mut_subtypes))
}

mcols(vcf_mini)$mut_type <- get_types(vcf_mini)
mcols(vcf_mini)$mut_subtype <-  get_subtypes(vcf_mini)
mcols(vcf)$mut_type <- get_types(vcf)
mcols(vcf)$mut_subtype <-  get_subtypes(vcf)

# distribucia mutacii na chromozomoch
get_distribution <- function(vcf_data, chromosome="all", subtypes=FALSE){
  chromosomes <- seqnames(rowRanges(vcf_data))
  if(!subtypes) types <- mcols(vcf_data)$mut_type
  else types <- mcols(vcf_data)$mut_subtype
  df <- data.frame(chromosomes, types)
  if (chromosome != "all") df %<>% filter(chromosomes == chromosome)
  ggplot(df, aes(x = chromosomes, fill = types)) +
    geom_bar(position = "stack", color = "black") +
    labs(title = "Variant Distribution on Chromosomes",
         x = "Chromosome",
         y = "Variant count") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

get_distribution(vcf, chromosome="chr1")
get_distribution(vcf, subtypes=TRUE)

#alelicka frekvencia 
get_frequency <- function(vcf_data){
  ad <- geno(vcf_data)$AD
  dp <- geno(vcf_data)$DP
  alt_AF <- map2_dbl(ad, dp, ~ {
    if (!is.null(.x) && !is.null(.y) && .y > 0) {
      1 - .x[1] / .y
    } else {NA}
  })
  results <- data.frame(Mutation =seq_along(alt_AF), alt_AF)
  ggplot(results) +
    geom_col(aes(x = Mutation, y = alt_AF), fill = "skyblue") +
    labs(title = "Allele Frequency",
         x = "Mutation",
         y = "Frequency") + theme_minimal()
}

get_frequency(vcf_mini)



