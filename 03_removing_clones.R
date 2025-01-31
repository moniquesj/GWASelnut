library(vcfR)
library(adegenet)
library(ape)
library(poppr)
library(plyr)
library(reshape2)


############################################################################
#### Getting a list of samples that are too similar and must be removed ####
############################################################################
setwd("/output/tassel_distance_matrix")

# Before importing the distance matrix done with Tassel, the first 5 rows must be deleted on Notepad
dist_matrix <- read.table("distmatrix_population.snps.hz.tombul.filtered.leaf.homozygous.txt",row.names = 1)

# Convert the data frame to a matrix
dist_matrix <- as.matrix(dist_matrix)
# Matching column names to the row names
colnames(dist_matrix) <- rownames(dist_matrix)
# Change the matrix from a wide format to a long format 
dist_matrix_melt <- melt(dist_matrix)


# Importing taxa summary
taxa_summary <- read.table("taxasummary_population.snps.hz.tombul.filtered.leaf.homozygous.txt", sep = "\t", header = T)

# Using mapvalues, do the following:
# Create a new column called Var_1_proportion_missing inside dist_matrix_melt
# Look inside dist_matrix_melt$Var1 
# Check if dist_matrix_melt$Var1 and taxa_summary$Taxa.name are the same
# If they are the same, paste the taxa_summary$Proportion.Missing in dist_matrix_melt$Var_1_proportion_missing

# Everything has to be character
dist_matrix_melt$Var1 <- as.character(dist_matrix_melt$Var1)
dist_matrix_melt$Var2 <- as.character(dist_matrix_melt$Var2)
taxa_summary$Taxa.Name <- as.character(taxa_summary$Taxa.Name)

dist_matrix_melt$Var_1_proportion_missing <- mapvalues(
  dist_matrix_melt$Var1,
  from = taxa_summary$Taxa.Name,
  to = as.character(taxa_summary$Proportion.Missing),
  warn_missing = TRUE
)

dist_matrix_melt$Var_2_proportion_missing <- mapvalues(
  dist_matrix_melt$Var2,
  from = taxa_summary$Taxa.Name,
  to = as.character(taxa_summary$Proportion.Missing),
  warn_missing = TRUE
)

#changing the name of melt_matrix to data
data<-dist_matrix_melt

# Iterate over each row
for (i in 1:nrow(data)) {
  column_one <- data[i, 1]
  column_two <- data[i, 2]
  column_three <- as.numeric(data[i, 3])
  column_four <- as.numeric(data[i, 4])
  column_five <- as.numeric(data[i, 5])
  
  if (column_three < 0.001 & column_four > column_five) {
    data[i, "sample_to_be_deleted"] <- column_one
  } else if (column_three < 0.001 & column_four < column_five) {
    data[i, "sample_to_be_deleted"] <- column_two
  } else {
    cat("Row:", paste(data[i, ], collapse = ", "), "\n")
    data[i, "sample_to_be_deleted"] <- "none"
  }
}

# Which samples need to be kept now?
# Extracting the column from the table
sample_to_be_deleted <- data$sample_to_be_deleted

# Keep only unique values
unique_values_0.001 <- unique(sample_to_be_deleted)

# Remove the value 'none'
unique_values_0.001 <- unique_values_0.001[unique_values_0.001 != 'none']

# List of samples to be deleted from the vcf file on Tassel
# Write only selected_data$Taxa.Name
write.table(unique_values_0.001, file = "samples_to_be_deleted_0.001.txt", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

