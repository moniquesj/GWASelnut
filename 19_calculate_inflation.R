library(readr)

setwd("D:/OneDrive - Scuola Superiore Sant'Anna/BioLabs/Projects/Hazelnut/Data analysis/Second_sequencing/output/gwas_seed_traits/PCA.total = 3")

# Import the output of GWAS for nut average size
gwas_nut_size <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.average.size.csv", 
                               col_types = cols(Chr = col_number(), 
                               Pos = col_number(), P.value = col_number(), 
                               MAF = col_number(), nobs = col_number(), 
                               `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_size <- gwas_nut_size$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_size <- qchisq(1-p_values_nut_size,1)

# Calculate lambda gc (λgc)
median(chisq_nut_size)/qchisq(0.5,1) # 0.8551531

# Import the output of GWAS for nut average size
gwas_nut_size <- read_csv("GAPIT.Association.GWAS_Results.BLINK.Nut.average.size.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_size <- gwas_nut_size$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_size <- qchisq(1-p_values_nut_size,1)

# Calculate lambda gc (λgc)
median(chisq_nut_size)/qchisq(0.5,1) # 0.8868401

# Import the output of GWAS for nut average size
gwas_nut_size <- read_csv("GAPIT.Association.GWAS_Results.GLM.Nut.average.size.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_size <- gwas_nut_size$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_size <- qchisq(1-p_values_nut_size,1)

# Calculate lambda gc (λgc)
median(chisq_nut_size)/qchisq(0.5,1) # 2.141191

# Import the output of GWAS for nut average size
gwas_nut_size <- read_csv("GAPIT.Association.GWAS_Results.MLM.Nut.average.size.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_size <- gwas_nut_size$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_size <- qchisq(1-p_values_nut_size,1)

# Calculate lambda gc (λgc)
median(chisq_nut_size)/qchisq(0.5,1) # 0.9660003

######################################################################################################
# Import the output of GWAS for nut perimeter
gwas_nut_perimeter <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.perimeter.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_perimeter <- gwas_nut_perimeter$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_perimeter <- qchisq(1-p_values_nut_perimeter,1)

# Calculate lambda gc (λgc)
median(chisq_nut_perimeter)/qchisq(0.5,1) # 0.7338432

gwas_nut_perimeter <- read_csv("GAPIT.Association.GWAS_Results.BLINK.Nut.perimeter.csv", 
                          col_types = cols(Chr = col_number(), 
                                           Pos = col_number(), P.value = col_number(), 
                                           MAF = col_number(), nobs = col_number(), 
                                           `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_perimeter <- gwas_nut_perimeter$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_perimeter <- qchisq(1-p_values_nut_perimeter,1)

# Calculate lambda gc (λgc)
median(chisq_nut_perimeter)/qchisq(0.5,1) # 1.040195

gwas_nut_perimeter <- read_csv("GAPIT.Association.GWAS_Results.GLM.Nut.perimeter.csv", 
                               col_types = cols(Chr = col_number(), 
                                                Pos = col_number(), P.value = col_number(), 
                                                MAF = col_number(), nobs = col_number(), 
                                                `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_perimeter <- gwas_nut_perimeter$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_perimeter <- qchisq(1-p_values_nut_perimeter,1)

# Calculate lambda gc (λgc)
median(chisq_nut_perimeter)/qchisq(0.5,1) #  1.983642


gwas_nut_perimeter <- read_csv("GAPIT.Association.GWAS_Results.MLM.Nut.perimeter.csv", 
                               col_types = cols(Chr = col_number(), 
                                                Pos = col_number(), P.value = col_number(), 
                                                MAF = col_number(), nobs = col_number(), 
                                                `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_perimeter <- gwas_nut_perimeter$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_perimeter <- qchisq(1-p_values_nut_perimeter,1)

# Calculate lambda gc (λgc)
median(chisq_nut_perimeter)/qchisq(0.5,1) #  1.042875

####################################################################################################

######################################################################################################
# Import the output of GWAS for nut maximum caliber
gwas_nut_maximum_caliber <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.maximum.caliber.csv", 
                               col_types = cols(Chr = col_number(), 
                                                Pos = col_number(), P.value = col_number(), 
                                                MAF = col_number(), nobs = col_number(), 
                                                `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_maximum_caliber <- gwas_nut_maximum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_maximum_caliber <- qchisq(1-p_values_nut_maximum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_maximum_caliber)/qchisq(0.5,1) # 0.8713541

gwas_nut_maximum_caliber <- read_csv("GAPIT.Association.GWAS_Results.BLINK.Nut.maximum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_maximum_caliber <- gwas_nut_maximum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_maximum_caliber <- qchisq(1-p_values_nut_maximum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_maximum_caliber)/qchisq(0.5,1) # 1.157148

gwas_nut_maximum_caliber <- read_csv("GAPIT.Association.GWAS_Results.GLM.Nut.maximum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_maximum_caliber <- gwas_nut_maximum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_maximum_caliber <- qchisq(1-p_values_nut_maximum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_maximum_caliber)/qchisq(0.5,1) # 2.509461

gwas_nut_maximum_caliber <- read_csv("GAPIT.Association.GWAS_Results.MLM.Nut.maximum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_maximum_caliber <- gwas_nut_maximum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_maximum_caliber <- qchisq(1-p_values_nut_maximum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_maximum_caliber)/qchisq(0.5,1) # 0.9825204


######################################################################################################
# Import the output of GWAS for nut minimum caliber
gwas_nut_minimum_caliber <- read_csv("GAPIT.Association.GWAS_Results.FarmCPU.Nut.minimum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_minimmum_caliber <- gwas_nut_minimum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_minimum_caliber <- qchisq(1-p_values_nut_minimmum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_minimum_caliber)/qchisq(0.5,1) # 0.7525167

gwas_nut_minimum_caliber <- read_csv("GAPIT.Association.GWAS_Results.BLINK.Nut.minimum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_minimmum_caliber <- gwas_nut_minimum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_minimum_caliber <- qchisq(1-p_values_nut_minimmum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_minimum_caliber)/qchisq(0.5,1) # 0.8282933

gwas_nut_minimum_caliber <- read_csv("GAPIT.Association.GWAS_Results.GLM.Nut.minimum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_minimmum_caliber <- gwas_nut_minimum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_minimum_caliber <- qchisq(1-p_values_nut_minimmum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_minimum_caliber)/qchisq(0.5,1) # 2.68128

gwas_nut_minimum_caliber <- read_csv("GAPIT.Association.GWAS_Results.MLM.Nut.minimum.caliber.csv", 
                                     col_types = cols(Chr = col_number(), 
                                                      Pos = col_number(), P.value = col_number(), 
                                                      MAF = col_number(), nobs = col_number(), 
                                                      `H&B.P.Value` = col_number(), Effect = col_number()))

p_values_nut_minimmum_caliber <- gwas_nut_minimum_caliber$P.value

# For p-values, calculate chi-squared statistic
chisq_nut_minimum_caliber <- qchisq(1-p_values_nut_minimmum_caliber,1)

# Calculate lambda gc (λgc)
median(chisq_nut_minimum_caliber)/qchisq(0.5,1) # 0.947052

