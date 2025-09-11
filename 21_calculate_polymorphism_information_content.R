# Calculating Polymorphism Information Content
#  For a biallelic SNP, the Polymorphism Information Content (PIC) is calculated as:
#  PIC= 1 - (p^2 + q^2) - 2 * (p^2 * q^2)
# where p is the major allele frequency and q is the minor allele frequency  
# This formula was developed by Botstein et al., 1980 in the paper Construction of a genetic linkage map in man using restriction fragment length polymorphisms. American Journal of Human Genetics, 32, 314–331.

library(readr)
library(dplyr)

####################################
########### Nut area PVE ###########
####################################
# Importing GWAS output
nut_area_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.average.size.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))
nut_area_gwas$Trait <- "Area"

# Keep only the SNPs that were significantly associated to the trait
significant_snps_area <- nut_area_gwas[order(nut_area_gwas$P.value), ][1:3, ]

# Calculate PIC
significant_snps_area$PIC <- with(significant_snps_area, {
  q <- MAF
  p <- 1 - q
  1 - (p^2 + q^2) - 2 * (p^2 * q^2)
})



####################################
########### Nut maximum caliber PVE 
####################################
# Importing GWAS output
nut_max_caliber_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.maximum.caliber.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))
nut_max_caliber_gwas$Trait <- "Maximum caliber"

# Keep only the SNPs that were significantly associated to the trait
significant_snps_max_caliber <- nut_max_caliber_gwas[order(nut_max_caliber_gwas$P.value), ][1:2, ]

# Calculate PIC
significant_snps_max_caliber$PIC <- with(significant_snps_max_caliber, {
  q <- MAF
  p <- 1 - q
  1 - (p^2 + q^2) - 2 * (p^2 * q^2)
})

####################################
########### Nut minimum caliber PVE 
####################################
nut_min_caliber_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.minimum.caliber.csv", 
                                 col_types = cols(Chr = col_number(), 
                                                  Pos = col_number(), P.value = col_number(), 
                                                  MAF = col_number(), nobs = col_number(), 
                                                  `H&B.P.Value` = col_number(), Effect = col_number()))
nut_min_caliber_gwas$Trait <- "Minimum caliber"

# Keep only the SNPs that were significantly associated to the trait
significant_snps_min_caliber <- nut_min_caliber_gwas[order(nut_min_caliber_gwas$P.value), ][1:3, ]

# Calculate PIC
significant_snps_min_caliber$PIC <- with(significant_snps_min_caliber, {
  q <- MAF
  p <- 1 - q
  1 - (p^2 + q^2) - 2 * (p^2 * q^2)
})

####################################
########### Nut perimeter PVE 
####################################
nut_perimeter_gwas <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.perimeter.csv", 
                                 col_types = cols(Chr = col_number(), 
                                                  Pos = col_number(), P.value = col_number(), 
                                                  MAF = col_number(), nobs = col_number(), 
                                                  `H&B.P.Value` = col_number(), Effect = col_number()))

nut_perimeter_gwas$Trait <- "Perimeter"

# Keep only the SNPs that were significantly associated to the trait
significant_snps_perimeter <- nut_perimeter_gwas[order(nut_perimeter_gwas$P.value), ][1:2, ]

significant_snps_perimeter$PIC <- with(significant_snps_perimeter, {
  q <- MAF
  p <- 1 - q
  1 - (p^2 + q^2) - 2 * (p^2 * q^2)
})

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
write.csv(all_significant_snps, "PIC_all_significant_snps.csv", row.names = FALSE)
