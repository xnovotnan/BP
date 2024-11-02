library(ggplot2)
library(BiocManager)
library(VariantAnnotation)
library(tidyr)

vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new.vcf"
vcf_file_mini <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N_new_mini.vcf"
vcf <- readVcf(vcf_file, "t2t")
vcf_mini <- readVcf(vcf_file_mini, "t2t")


create_dataframe <- function(vcf){
  chromosomes <-seqnames(rowRanges(vcf))
  positions <- start(rowRanges(vcf))
  ref <- rowRanges(vcf)$REF
  alt <- rowRanges(vcf)$ALT
  alt <- CharacterList(alt)
  alt <- unstrsplit(alt, sep = ",")
  type <- character(length(ref))
  subtype <- character(length(ref))
  
  vcf_data <- data.frame(
    Chromosome = as.character(chromosomes),
    Position = positions,
    Ref = as.character(ref),
    Alt = alt,
    stringsAsFactors = FALSE, 
    Type = type,
    Subtype = subtype
  )
  return(vcf_data)
}

vcf_data <- create_dataframe(vcf_mini)


# Alelicke frekvencie ref/alt + graf
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
    geom_col(aes(x = Mutation, y = ref_AF, fill = "Reference Allele Frequency")) +
    geom_col(aes(x = Mutation, y = alt_AF, fill = "Alternative Allele Frequency")) +
    labs(title = "Reference and Alternative Allele Frequencies",
         x = "Mutation",
         y = "Frequency") 
  }

plot <- get_AF(vcf_mini)
plot
ggsave("allele_frequencies.png", plot = plot, width = 10, height = 6)



# Distribucia mutacii na chromozomoch
get_chrom_hist <- function(vcf, chromosome="all"){
  if (chromosome=="all"){
    chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    vcf_data$Chromosome <- factor(vcf_data$Chromosome, levels = chromosome_order)
    ggplot(data = vcf_data, aes(x = Chromosome)) +
      geom_bar(fill="skyblue", color = "black") + 
      labs(title = "Variant Distribution on Chromosomes",
           x = "Chromosome",
           y = "Variant count") + theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  else{
    chr <- subset(vcf_data, Chromosome == chromosome)
    chr$Segment <- cut(chr$Position, breaks = seq(0, max(chr$Position), by = 10000000), labels = FALSE)
    
    ggplot(chr, aes(x = Segment)) + 
      geom_bar(fill = "skyblue", color = "black") +
      labs(title = paste("Variant Distribution on Chromosome", chromosome),
           x = "Chromosome Segment (10 million base pairs per bin)",
           y = "Variant count") + theme_minimal()
  }
}

get_chrom_hist(vcf_mini, chromosome = "chr1")
get_chrom_hist(vcf_mini)


# Klasifikovanie mutacii 
classify_mutations <- function(vcf_data){
  purines <- c("A", "G")
  pyrimidines <- c("C", "T")
  
  #vcf_data vytvorene cez create_dataframe
  for (i in seq_along(vcf_data$Ref)){
    cur_ref = vcf_data$Ref[i]
    cur_alt = vcf_data$Alt[i]
    
    if(nchar(cur_ref) == 1 && nchar(cur_alt) == 1){
      vcf_data$Type[i] <- "SNP"
      if ((cur_ref %in% purines && cur_alt %in% purines) || 
          (cur_ref %in% pyrimidines && cur_alt %in% pyrimidines)) {
        vcf_data$Subtype[i] <- "Transition" 
      }else {
        vcf_data$Subtype[i] <- "Transversion"
      }
    }else if(nchar(cur_ref) > nchar(cur_alt)){
      vcf_data$Type[i] <- "INDEL"
      vcf_data$Subtype[i] <- "Deletion"
    }else if(nchar(cur_ref) < nchar(cur_alt)){
      vcf_data$Type[i] <- "INDEL"
      vcf_data$Subtype[i] <- "Insertion"
    }else {
      vcf_data$Type[i] <- "other"
    }
  }
  
  return(vcf_data)
}

vcf_data <- classify_mutations(vcf_data)
vcf_data['Type']




# da sa v R ulozit datova struktura ako singleton? alebo ulozit globalne a robit kopie? 

# Distribucia mutacii na chromozomoch s typmi mutacii
get_chrom_hist2 <- function(vcf, chromosome="all", include_types=FALSE){
  chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  
  if (chromosome == "all") {
    vcf_data$Chromosome <- factor(vcf_data$Chromosome, levels = chromosome_order)
    
    if (include_types) {
      ggplot(data = vcf_data, aes(x = Chromosome, fill = Type)) +
        geom_bar(color = "black") + 
        labs(title = "Variant Distribution on Chromosomes",
             x = "Chromosome",
             y = "Variant count") + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_fill_brewer(palette = "Set3") 
    } else {
      ggplot(data = vcf_data, aes(x = Chromosome)) +
        geom_bar(fill="skyblue", color = "black") + 
        labs(title = "Variant Distribution on Chromosomes",
             x = "Chromosome",
             y = "Variant count") + 
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    }
  } else {
    chr <- subset(vcf_data, Chromosome == chromosome)
    chr$Segment <- cut(chr$Position, breaks = seq(0, max(chr$Position), by = 10000000), labels = FALSE)
    
    if (include_types) {
      ggplot(chr, aes(x = Segment, fill = Type)) + 
        geom_bar(color = "black") +
        labs(title = paste("Variant Distribution on Chromosome", chromosome),
             x = "Chromosome Segment (10 million base pairs per bin)",
             y = "Variant count") + 
        theme_minimal() +
        scale_fill_brewer(palette = "Set3")
    } else {
      ggplot(chr, aes(x = Segment)) + 
        geom_bar(fill = "skyblue", color = "black") +
        labs(title = paste("Variant Distribution on Chromosome", chromosome),
             x = "Chromosome Segment (10 million base pairs per bin)",
             y = "Variant count") + 
        theme_minimal()
    }
  }
}
get_chrom_hist2(vcf_mini, include_types=TRUE)
get_chrom_hist2(vcf_mini, chromosome = "chr1", include_types=TRUE)
