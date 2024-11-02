library(ggplot2)
library(BiocManager)
library(VariantAnnotation)
library(tidyr)

vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new.vcf"
vcf_file_mini <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new_mini.vcf"
vcf <- readVcf(vcf_file, "t2t")
vcf_mini <- readVcf(vcf_file_mini, "t2t")

get_AF <- function(vcf){
  ad <- geno(vcf)$AD
  dp <- geno(vcf)$DP
  ref_AF <- numeric(length(ad))
  alt_AF <- numeric(length(ad))
  
  for (i in seq_along(ad)) {
    if (!is.null(ad[[i]]) && !is.null(dp[[i]])) {
      if (dp[[i]] > 0) { 
        af <- ad[[i]][1] / dp[[i]]
        ref_AF[i] <- af
        alt_AF[i] <- 1 - af
      } else {
        ref_AF[i] <- NA # 0 či NA
        alt_AF[i] <- NA
      }
    } else {
      ref_AF[i] <- NA
      alt_AF[i] <- NA
    }
  }
  
  results <- data.frame(
    Mutation = rownames(vcf),
    ref_AF = ref_AF,
    alt_AF = alt_AF
  )

  ggplot(results) +
    geom_col(aes(x = Mutation, y = ref_AF, fill = "Blue")) +
    geom_col(aes(x = Mutation, y = alt_AF, fill = "Red")) +
    labs(title = "Reference and Alternative Allele Frequencies",
         x = "Mutation",
         y = "Frequency") 
  }

get_AF(vcf_mini)



