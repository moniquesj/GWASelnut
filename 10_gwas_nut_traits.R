library(readr)
library(GAPIT)
library(PerformanceAnalytics)

# Set wd
setwd("/datasets/gwas")

# Importing the phenotype file
pheno <- read_delim("filtered.data.nut.traits.gwas.txt",
                    delim = "\t", escape_double = FALSE,
                    trim_ws = TRUE)
pheno <- as.data.frame(pheno)
pheno[, 5:11] <- lapply(pheno[, 5:11], as.numeric)

# Keeping only nut traits and Sample
columns_to_keep <- c("Sample", "Nut.average.size", "Nut.perimeter", 
                     "Nut.maximum.caliber", "Nut.minimum.caliber")

# Subset the data frame to keep only the specified columns
pheno_nut <- pheno[, columns_to_keep, drop = FALSE]

# Adjust first column name
colnames(pheno_nut)[1]<-"Taxa"

######################################
########### GWAS #####################
######################################

#Importing the HapMap diploid file with 151 samples (the vcf file was converted using Tassel)
genotypes <- read.delim("population.snps.hz.tombul.filtered.leaf.nut.gwas.hmp.txt", header = F)

#PCA=3
setwd("gwas_nut_traits/PCA.total = 3")
myGAPIT=GAPIT(
  Y=pheno_nut, #first column is ID
  G=genotypes,
  PCA.total=3, 
  model=c("GLM", "MLM", "Blink","FarmCPU"),
  Multiple_analysis=TRUE,
  Random.model = FALSE)

