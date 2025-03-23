library(tidyverse)
library(magrittr)
library(patchwork)
library(ggridges)

# -----------------------------------
# PREDSPRACOVANIE DAT A ZAKLADNA VIZUALIZACIA VCF SUBORU
# -----------------------------------

# Funkcia prepare_data() sluzi na spracovanie VCF suboru - vytvorenie tibble, 
# vymazanie hlavicky, vymazanie pozorovani s qual < 200 a s GT 0/0, 
# klasifikovanie typu a subtypu, vypocet alelickej frekvencie

prepare_data <- function(vcf_file){
  suppressWarnings(suppressMessages({
    sink("/dev/null")
    incProgress(0.1, detail = "Reading file...")
    header_line <- grep("^#CHROM", readLines(vcf_file)) 
    vcf_tibble <- read_delim(vcf_file, delim = "\t", skip = header_line - 1)
    incProgress(0.3, detail = "Processing data...")
    vcf_tibble %<>% 
      rename("VALUES" = last(colnames(.)), "CHROM" = first(colnames(.))) %>%
      mutate(
        GT = str_split(VALUES, ":") %>% map_chr(~ .x[1]),
        AD = str_split(VALUES, ":") %>% map_chr(~ .x[2]) %>% str_split(",") %>% map_int(~ as.integer(.x[1])),
        DP = str_split(VALUES, ":") %>% map_int(~ as.integer(.x[3])),
        TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"SNV", "INDEL")) %>%
      filter(QUAL > 200, 
             GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"),
             DP>0,
             CHROM != "chrM")
    incProgress(0.7, detail = "Categorizing mutations...")
    purines <- c("A", "G")
    pyrimidines <- c("C", "T")
    vcf_tibble %<>% mutate(
      SUBTYPE = ifelse(TYPE == "SNV", 
                       ifelse((REF %in% purines & ALT %in% purines) | (REF %in% pyrimidines & ALT %in% pyrimidines), "Transition", "Transversion"), 
                       ifelse(nchar(REF) > nchar(ALT),"Deletion", "Insertion")),
      AF = round(1 - AD/DP, 2)
    )
    incProgress(0.9, detail = "Completing preprocessing...")
    vcf_tibble$CHROM <- factor(vcf_tibble$CHROM, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
    vcf_tibble %<>% select(c("CHROM", "POS", "REF", "ALT", "TYPE", "SUBTYPE", "QUAL", "GT", "AD", "DP", "AF"))
    sink(NULL)
    vcf_tibble
  }))
}

# Funkcia plot_summary() sluzi na vytvorenie zakladnych vizualizacii 
# distribucie hodnot jednotlivych atributov - kvality, hlbky citania,
# alelickej frekvencie, alelickej hlbky, pocet mutacii na jednotlivych 
# chromozomoch, pocty jednotlivych typov a subtypov mutacii
plot_summary <- function(data){
  boxplot_qual <- ggplot(data, aes(x = "", y = QUAL)) + 
    geom_boxplot(fill = "skyblue", alpha = 0.5) + 
    labs(title = "Distribution of Quality", x = "", y="") 
  
  boxplot_dp <- ggplot(data, aes(x = "", y = DP)) + 
    geom_boxplot(fill = "lightgreen", alpha = 0.5) + 
    labs(title = "Distribution of Read Depth", x = "", y="") 
  
  boxplot_af <- ggplot(data, aes(x = "", y = AF)) + 
    geom_boxplot(fill = "salmon", alpha = 0.5) + 
    labs(title = "Distribution of Allelic Frequency", x = "", y="")
  
  boxplot_ad <- ggplot(data, aes(x = "", y = AD)) + 
    geom_boxplot(fill = "orange", alpha = 0.5) + 
    labs(title = "Distribution of Allelic Depth", x = "", y="") 
  
  hist_chrom <- ggplot(data, aes(x = CHROM)) +
    geom_bar(fill = "purple", alpha = 0.5) +
    labs(title = "Chromosomes Counts", x = "", y="") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
  
  hist_type <- ggplot(data, aes(x = TYPE)) +
    geom_bar(fill = "salmon", alpha = 0.5) + 
    labs(title = "Type Counts", x = "", y="") 
  
  hist_subtype <- ggplot(data, aes(x = SUBTYPE)) +
    geom_bar(fill = "skyblue", alpha = 0.5) + 
    labs(title = "Subtype Counts", x = "", y="") 
  
  (hist_chrom | hist_type | hist_subtype) / 
    (boxplot_qual | boxplot_dp | boxplot_af | boxplot_ad)
    
}


# plot_summary <- function(data){
#   boxplot_qual <- plot_ly(data, y = ~QUAL, type = "box", boxmean = "sd", marker = list(color = 'skyblue')) %>%
#     layout(title = "Boxplot of Quality", yaxis = list(title = ""), xaxis = list(title = ""))
#   
#   boxplot_dp <- plot_ly(data, y = ~DP, type = "box", boxmean = "sd", marker = list(color = 'lightgreen')) %>%
#     layout(title = "Boxplot of Read Depth", yaxis = list(title = ""), xaxis = list(title = ""))
#   
#   boxplot_af <- plot_ly(data, y = ~AF, type = "box", boxmean = "sd", marker = list(color = 'salmon')) %>%
#     layout(title = "Boxplot of Allelic Frequency", yaxis = list(title = ""), xaxis = list(title = ""))
#   
#   boxplot_ad <- plot_ly(data, y = ~AD, type = "box", boxmean = "sd", marker = list(color = 'orange')) %>%
#     layout(title = "Boxplot of Allelic Depth", yaxis = list(title = ""), xaxis = list(title = ""))
#   
#   hist_chrom <- plot_ly(data, x = ~CHROM, type = "histogram", marker = list(color = 'purple', opacity = 0.5)) %>%
#     layout(title = "Chromosome Counts", xaxis = list(title = ""), yaxis = list(title = ""))
#   
#   hist_type <- plot_ly(data, x = ~TYPE, type = "histogram", marker = list(color = 'pink', opacity = 0.5)) %>%
#     layout(title = "Type Counts", xaxis = list(title = ""), yaxis = list(title = ""))
#   
#   hist_subtype <- plot_ly(data, x = ~SUBTYPE, type = "histogram", marker = list(color = 'lightblue', opacity = 0.5)) %>%
#     layout(title = "Subtype Counts", xaxis = list(title = ""), yaxis = list(title = ""))
#   
#   subplot(hist_chrom, hist_type, hist_subtype, nrows = 1) %>%
#     subplot(boxplot_qual, boxplot_dp, boxplot_af, boxplot_ad, nrows = 2)
# }


# prepare_data <- function(vcf_file){
#   header_line <- grep("^#CHROM", readLines(vcf_file)) 
#   vcf_tibble <- read_delim(vcf_file, delim = "\t", skip = header_line - 1)
#   vcf_tibble %<>% 
#     rename("VALUES" = last(colnames(.)), "CHROM" = first(colnames(.))) %>%
#     mutate(
#       GT = str_split(VALUES, ":") %>% map_chr(~ .x[1]),
#       AD = str_split(VALUES, ":") %>% map_chr(~ .x[2]) %>% str_split(",") %>% map_int(~ as.integer(.x[1])),
#       DP = str_split(VALUES, ":") %>% map_int(~ as.integer(.x[3])),
#       TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"SNV", "INDEL")) %>%
#     filter(QUAL > 200, 
#            GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"),
#            DP>0,
#            CHROM != "chrM")
#   purines <- c("A", "G")
#   pyrimidines <- c("C", "T")
#   vcf_tibble %<>% mutate(
#     SUBTYPE = ifelse(TYPE == "SNV", 
#                      ifelse((REF %in% purines & ALT %in% purines) | (REF %in% pyrimidines & ALT %in% pyrimidines), "Transition", "Transversion"), 
#                      ifelse(nchar(REF) > nchar(ALT),"Deletion", "Insertion")),
#     AF = round(1 - AD/DP, 2)
#   )
#   vcf_tibble$CHROM <- factor(vcf_tibble$CHROM, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
#   vcf_tibble %<>% select(c("CHROM", "POS", "REF", "ALT", "TYPE", "SUBTYPE", "QUAL", "GT", "AD", "DP", "AF"))
#   vcf_tibble
# }
# a <- prepare_data("/Users/macbook/Documents/BP/data/Lynch.2526.01.N.vcf")
