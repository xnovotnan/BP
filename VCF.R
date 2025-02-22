library(VariantAnnotation)
library(BiocManager)
library(tidyverse)
library(magrittr)
library(patchwork)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(plotly)
library(circlize)
library(ComplexHeatmap)

setwd("/Users/macbook/Documents/BP")
vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N.vcf"

theme_set(theme_minimal() + 
            theme(
              plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
              legend.position = "bottom"
            ))
options(scipen = 999)


# -----------------------------------
# PREDSPRACOVANIE DAT
# -----------------------------------
# vytvorenie tibble, vymazanie hlavicky a mutacii s qual < 200 a s GT 0/0, 
# klasifikovanie typu a subtypu, vypocet alelickej frekvencie
prepare_data <- function(vcf_file){
  header_line <- grep("^#CHROM", suppressMessages(readLines(vcf_file))) 
  vcf_tibble <- suppressMessages(read_delim(vcf_file, delim = "\t", skip = header_line - 1))
  vcf_tibble %<>% 
    rename("VALUES" = last(colnames(.)), "CHROM" = first(colnames(.))) %>%
    mutate(
      GT = str_split(VALUES, ":") %>% map_chr(~ .x[1]),
      AD = str_split(VALUES, ":") %>% map_chr(~ .x[2]) %>% str_split(",") %>% map_int(~ as.integer(.x[1])),
      DP = str_split(VALUES, ":") %>% map_int(~ as.integer(.x[3])),
      TYPE = ifelse(nchar(REF) == 1 & nchar(ALT) == 1,"SNP", "INDEL")) %>%
    filter(QUAL > 200, 
           GT %in% c("1|1", "0|1", "1|0", "1/1", "0/1", "1/0"),
           DP>0)

  purines <- c("A", "G")
  pyrimidines <- c("C", "T")
  vcf_tibble %<>% mutate(
    SUBTYPE = ifelse(TYPE == "SNP", 
                     ifelse((REF %in% purines & ALT %in% purines) | (REF %in% pyrimidines & ALT %in% pyrimidines), "Transition", "Transversion"), 
                     ifelse(nchar(REF) > nchar(ALT),"Deletion", "Insertion")),
    AF = 1 - AD/DP
  )
  
  vcf_tibble$CHROM <- factor(vcf_tibble$CHROM, levels = c(paste0("chr", 1:22), "chrX", "chrY", "chrM"))
  vcf_tibble %<>% select(c("CHROM", "POS", "REF", "ALT", "TYPE", "SUBTYPE", "QUAL", "GT", "AD", "DP", "AF"))
  vcf_tibble
}

vcf_complete <- prepare_data(vcf_file)
vcf_complete

# -----------------------------------
# GRAFY KVALITY A HLBKY CITANI
# -----------------------------------

# Graf hodnôt kvality
plot_QUAL <- function(vcf_tibble) {
  p <- ggplot(vcf_tibble, aes(x=seq_along(QUAL), y=QUAL)) +
    geom_point() + 
    labs(title = "Quality along the VCF file", 
         x = "Position", 
         y="Quality") +
    theme_minimal()
  ggplotly(p)
}
plot_QUAL(vcf_complete)


# Violin plot kvality na jednotlivých chromozómoch
qual_on_chroms <- function(vcf_tibble, include_M=FALSE){
  if(!include_M){
    vcf_tibble %<>% filter(!CHROM %in% c("chrM"))
  }
  p <- ggplot(vcf_tibble, aes(x=CHROM, y=QUAL, fill=CHROM)) +
    geom_violin() +
    labs(title = "Quality across chromosomes",
         x = "Chromosomes",
         y = "Quality") +
    theme_minimal() +
    theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggplotly(p)
}
qual_on_chroms(vcf_complete)
qual_on_chroms(vcf_complete, include_M = TRUE)


# Graf hĺbky čítaní 
read_depth <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x=DP)) +
    geom_bar(fill = "skyblue") + 
    labs(title = "Distribution of read depth",
         x = "",
         y = "Read depth") +
    theme_minimal() 
  ggplotly(p)
}
read_depth(vcf_complete)


# Read depth across chromosomes - violin plot
rd_on_chroms <- function(vcf_tibble, include_M=FALSE){
  if(!include_M){
    vcf_tibble %<>% filter(!CHROM %in% c("chrM"))
  }
  p <- ggplot(vcf_tibble, aes(x=CHROM, y=DP, fill=CHROM)) +
    geom_violin() +
    labs(title = "Read depth across chromosomes",
         x = "Chromosomes",
         y = "Read depth") +
    theme_minimal() +
    theme(legend.position="none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggplotly(p)
}
rd_on_chroms(vcf_complete)
rd_on_chroms(vcf_complete, include_M = TRUE)


# Hexbin QUAL vs Read depth
qual_vs_rd <- function(vcf_tibble){
  p <- ggplot(vcf_tibble, aes(x = DP, y = QUAL)) +
    geom_hex() + 
    #scale_x_log10() +
    theme_minimal()
  ggplotly(p)
}
qual_vs_rd(vcf_complete)


# -----------------------------------
# GRAFY ALELICKEJ FREKVENCIE
# -----------------------------------

# Graf alelickej frekvencie - denisty plot + hexbin
plot_frequency <- function(vcf_tibble) {
  p <- ggplot(vcf_tibble, aes(x=AF)) +
    geom_density(fill = "skyblue", alpha=0.5) +
    labs(title = "Allele Frequency",
         x = "Mutation",
         y = "Frequency") +
    theme_minimal()
  ggplotly(p)
}
plot_frequency(vcf_complete)

plot_frequency_hexbin <- function(vcf_tibble) {
  p <- ggplot(vcf_tibble, aes(x=seq_along(AF), y = AF)) +
    geom_hex() +
    labs(title = "Allele Frequency",
         x = "",
         y = "Frequency") +
    theme_minimal()
  ggplotly(p)
}
plot_frequency_hexbin(vcf_complete)


#Alelická frekvencia per contig #plotly protestuje
AF_per_contig <- function(vcf_tibble){
  ggplot(vcf_tibble, aes(x=AF, y=CHROM, fill=CHROM)) +
    geom_density_ridges() +
    labs(title = "Allele Frequency") +
    theme_minimal() +
    theme(legend.position="none")
}

AF_per_contig(vcf_complete)


#Alellická frekvencia per contig - výber contigov
AF_per_certain_contigs <- function(vcf_tibble, contigs){
  vcf_tibble %<>% filter(CHROM %in% contigs)
  
  p <- ggplot(vcf_tibble, aes(x=AF, color=CHROM, fill=CHROM)) +
    geom_density(alpha=0.5) +
    labs(title = "Allele Frequency",
         x = "Mutation",
         y = "Frequency") +
    theme_minimal() +
    theme(legend.position="none")
  ggplotly(p)
}
AF_per_certain_contigs(vcf_complete, c("chr1", "chrX"))


# -----------------------------------
# TYPY MUTACII A ICH DISTRIBUCIA
# -----------------------------------

# Distribúcia mutácií - barplot a pie chart (skoro)
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
    labs(y=NULL, x=NULL)
  
  
  p2 <- ggplot(vcf_tibble, aes(x=TYPE, fill = TYPE)) +
    geom_bar(color = "black") + 
    theme_minimal()+
    coord_polar()+
    theme(legend.position="none", 
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    labs(y=NULL,x=NULL)
  
  (p1 + p2) + 
    plot_layout(ncol = 2) + 
    plot_annotation(title = "Mutation Type Distribution") &
    theme(plot.title = element_text(hjust = 0.5))
}

mut_summary(vcf_complete)
mut_summary(vcf_complete, subtypes=TRUE)


# Distribúcia mutácií na jednotlivých chromozómoch 
# možnosť zvoliť úroveň klasifikácie mutácií aj konkrétny chromozóm
get_distribution <- function(vcf_tibble, chromosome="all", subtypes=FALSE){
  if (chromosome != "all") vcf_tibble %<>% filter(CHROM == chromosome) 
  if (subtypes) vcf_tibble$TYPE <- vcf_tibble$SUBTYPE
  
  p <- ggplot(vcf_tibble, 
         aes(x = CHROM, fill = TYPE)) +
    geom_bar(position = "stack", color = "black") +
    labs(title = "Variant Distribution on Chromosomes",
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)) +
    coord_flip()
  ggplotly(p)
}

get_distribution(vcf_complete, chromosome="chr1")
get_distribution(vcf_complete, chromosome="chr1", subtypes = TRUE)
get_distribution(vcf_complete, subtypes=TRUE)
get_distribution(vcf_complete, subtypes=FALSE)


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

mutation_summary_circle(vcf_complete)
mutation_summary_circle(vcf_complete, subtypes = TRUE)


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

lollipop(vcf_complete)

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

lollipop_subtypes(vcf_complete)

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

SNV_types(vcf_complete)


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
  
  ggplotly(p)
} 

mutation_heatmap(vcf_complete)


# -----------------------------------
# CIRCOS GRAFY
# -----------------------------------


#scattrrplot AF
circos.clear()
circos.par(cell.padding = c(0.02, 0, 0.02, 0))

circos.initialize(vcf_complete$CHROM, vcf_complete$POS)
circos.track(ylim = c(0, 1))
circos.trackPoints(vcf_complete$CHROM, vcf_complete$POS, vcf_complete$AF)



# heatmap AF
col_fun1 = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
circos.clear()
circos.heatmap(vcf_complete$AF, split= vcf_complete$CHROM, col = col_fun1)







runExample("01_hello")      # a histogram
runExample("02_text")       # tables and data frames
runExample("03_reactivity") # a reactive expression
runExample("04_mpg", display.mode = "showcase")        # global variables
runExample("05_sliders")    # slider bars
runExample("06_tabsets")    # tabbed panels
runExample("07_widgets")    # help text and submit buttons
runExample("10_download")   # file download wizard


summary(vcf_complete)
