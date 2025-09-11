# Calculating phenotypic variance explained (PVE) per SNP per trait

# The proportion of variance explained (PVE) by each association (SNP and seed trait) was calculated using the
# expression: VQTL = Vpheno, where VQTL = 2freq (1 – freq)effect^ 2 and Vpheno is the variance of the adjusted means for each trait.
# This formula comes from the paper Genome-wide association studies dissect the genetic architecture of seed shape and size in common bean 
# from Giordani et al. 2022

library(readr)
library(dplyr)


#importing the phenotype file
pheno <- read_delim("filtered.data.seed.traits.gwas.txt",
                    delim = "\t", escape_double = FALSE,
                    trim_ws = TRUE)
pheno <- as.data.frame(pheno)
pheno[, 5:11] <- lapply(pheno[, 5:11], as.numeric)

#keeping only Seed traits and VCF numbers
columns_to_keep <- c("All.VCF.numbers", "Nut.average.size", "Nut.perimeter", 
                     "Nut.maximum.caliber", "Nut.minimum.caliber")

# Subset the data frame to keep only the specified columns
pheno_seed <- pheno[, columns_to_keep, drop = FALSE]

#importing the list of the 151 samples effectively used for the GWAS
samples <- read_delim("samples_seed_traits_gwas.txt", 
                           delim = "\t", escape_double = FALSE, 
                           col_names = FALSE, trim_ws = TRUE)
#Keep in pheno_seed only the samples used for the GWAS
colnames(samples) <- "All.VCF.numbers"
pheno_seed <- pheno_seed[pheno_seed$All.VCF.numbers %in% samples$All.VCF.numbers, ]

####################################
########### Nut area PVE ###########
####################################
# Importing GWAS output
nut_area_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.average.size.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

# Keep only the SNPs that were significantly associated to the trait
significant_snps_area <- nut_area_gwas[order(nut_area_gwas$P.value), ][1:3, ]

# Calculate total variance of nut area trait
total_variance_area <- var(pheno_seed$Nut.average.size, na.rm = TRUE)
# Vqtl calculation 
significant_snps_area$Vqtl <- 2 * significant_snps_area$MAF * (1 - significant_snps_area$MAF) * (significant_snps_area$Effect ^ 2)
# Calculate PVE as percentage
significant_snps_area$PVE_percent <- (significant_snps_area$Vqtl / total_variance_area) * 100
significant_snps_area['Trait'] = 'Area'

####################################
########### Nut maximum caliber PVE 
####################################
# Importing GWAS output
nut_max_caliber_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.maximum.caliber.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

# Keep only the SNPs that were significantly associated to the trait
significant_snps_max_caliber <- nut_max_caliber_gwas[order(nut_max_caliber_gwas$P.value), ][1:2, ]

# Calculate total variance of maximum caliber trait
total_variance_max_caliber <- var(pheno_seed$Nut.maximum.caliber, na.rm = TRUE)

# Vqtl calculation 
significant_snps_max_caliber$Vqtl <- 2 * significant_snps_max_caliber$MAF * (1 - significant_snps_max_caliber$MAF) * (significant_snps_max_caliber$Effect ^ 2)
# Calculate PVE as percentage
significant_snps_max_caliber$PVE_percent <- (significant_snps_max_caliber$Vqtl / total_variance_max_caliber) * 100
significant_snps_max_caliber['Trait'] = 'Maximum caliber'

####################################
########### Nut minimum caliber PVE 
####################################
nut_min_caliber_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.minimum.caliber.csv", 
                                 col_types = cols(Chr = col_number(), 
                                                  Pos = col_number(), P.value = col_number(), 
                                                  MAF = col_number(), nobs = col_number(), 
                                                  `H&B.P.Value` = col_number(), Effect = col_number()))

# Keep only the SNPs that were significantly associated to the trait
significant_snps_min_caliber <- nut_min_caliber_gwas[order(nut_min_caliber_gwas$P.value), ][1:3, ]

# Calculate total variance of maximum caliber trait
total_variance_min_caliber <- var(pheno_seed$Nut.minimum.caliber, na.rm = TRUE)

# Vqtl calculation 
significant_snps_min_caliber$Vqtl <- 2 * significant_snps_min_caliber$MAF * (1 - significant_snps_min_caliber$MAF) * (significant_snps_min_caliber$Effect ^ 2)
# Calculate PVE as percentage
significant_snps_min_caliber$PVE_percent <- (significant_snps_min_caliber$Vqtl / total_variance_min_caliber) * 100
significant_snps_min_caliber['Trait'] = 'Minimum caliber'

####################################
########### Nut perimeter PVE 
####################################
nut_perimeter_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.perimeter.csv", 
                                 col_types = cols(Chr = col_number(), 
                                                  Pos = col_number(), P.value = col_number(), 
                                                  MAF = col_number(), nobs = col_number(), 
                                                  `H&B.P.Value` = col_number(), Effect = col_number()))

# Keep only the SNPs that were significantly associated to the trait
significant_snps_perimeter <- nut_perimeter_gwas[order(nut_perimeter_gwas$P.value), ][1:2, ]

# Calculate total variance of maximum caliber trait
total_variance_perimeter <- var(pheno_seed$Nut.perimeter, na.rm = TRUE)

# Vqtl calculation 
significant_snps_perimeter$Vqtl <- 2 * significant_snps_perimeter$MAF * (1 - significant_snps_perimeter$MAF) * (significant_snps_perimeter$Effect ^ 2)
# Calculate PVE as percentage
significant_snps_perimeter$PVE_percent <- (significant_snps_perimeter$Vqtl / total_variance_perimeter) * 100
significant_snps_perimeter['Trait'] = 'Perimeter'

#Create a table with all the results
all_significant_snps <- rbind(
  significant_snps_min_caliber,
  significant_snps_max_caliber,
  significant_snps_area,
  significant_snps_perimeter
)

# Make the column Trait the first of the df
all_significant_snps <- all_significant_snps %>%
  select(Trait, everything())

# Saving the results in a table
write.csv(all_significant_snps, "PVE_all_significant_snps.csv", row.names = FALSE)
