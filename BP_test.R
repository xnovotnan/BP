library(ggplot2)
library(BiocManager)
library(VariantAnnotation)

vcf_file <- "/Users/macbook/Documents/BP/data/Lynch.2526.01.N.vcf.gz"
vcf <- readVcf(vcf_file, "t2t")

chromosomes <-seqnames(rowRanges(vcf))
positions <- start(rowRanges(vcf))
types <- rowRanges(vcf)$REF  
alter <- rowRanges(vcf)$ALT
alter <- CharacterList(alter)
alter <- unstrsplit(alter, sep = ",")

vcf_data <- data.frame(
  Chromosome = as.character(chromosomes),
  Position = positions,
  Ref = as.character(types),
  Alt = alter,
  stringsAsFactors = FALSE
)

head(vcf_data)

chr1 <- subset(vcf_data, Chromosome == "chr1")
head(chr1)

ggplot(chr1, aes(x = Position)) + geom_density(fill = "skyblue") + 
  labs(title = "Density plot of variants on chr1", x = "Position", 
       y = "Density") + theme_minimal()

num_variants <- nrow(chr1)
max_pos <- tail(chr1,n=1)$Position

ggplot(chr1, aes(x = Position)) + 
  geom_histogram(binwidth = 1000000, fill = "skyblue", color = "black") +
  labs(title = "Histogram variantov na chromozóme chr1",
       x = "Genomic Position (bp)",
       y = "Počet variantov") + theme_minimal()

options(scipen = 999)
library(scales)
ggplot(chr1, aes(x = Position)) + 
  geom_density(fill = "lightblue") + 
  labs(title = "Density plot of variants on chr1", 
       x = "Position", 
       y = "Density") +
  xlim(0, max_pos) +
  scale_x_continuous(labels = comma) +
  theme_minimal()




  
  