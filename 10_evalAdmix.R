# plot evalAdmix results with K=5
setwd("/output/evalAdmix")
source("visFuns.R")

# Define a function to read evalAdmix files and perform the required operations
read_and_plot_evalAdmix <- function(K) {
  # Read fam, q, ped, and correlation matrix files
  fam <- read.table(paste0("population.snps.hz.tombul.filtered.pruned.leaf.for.plink.admixture.fam"), sep = "\t", header = FALSE)
  q <- read.table(paste0("population.snps.hz.tombul.filtered.pruned.leaf.for.plink.admixture.", K, ".Q"), stringsAsFactors = TRUE)
  ped <- read.table(paste0("population.snps.hz.tombul.filtered.pruned.leaf.for.plink.admixture.ped"), header = FALSE)
  r <- as.matrix(read.table(paste0("evaladmixOut.K", K, ".corres")))
  
  # Order fam file according to the sample names in the ped file
  ordered_fam <- fam[match(ped$V1, fam$V2), ]
  
  # Combine meta data and q data
  combined <- cbind(ordered_fam, q)
  
  # Order according to provenance categories
  combined_ordered <- combined[order(combined$V1), ]
  
  # Reorder the dataset, placing rows with provenance "Breeding" at the end
  combined_ordered_breeding <- rbind(
    combined_ordered[combined_ordered$V1 != "Breeding", ],
    combined_ordered[combined_ordered$V1 == "Breeding", ]
  )
  
  # Order according to population
  ord <- orderInds(pop = as.vector(combined_ordered_breeding[, 1]), q = q)
  
  # Open a PDF file for output with a unique name
  pdf(paste0("evaladmix_plot_K", K, ".pdf"))
  
  # Plot correlation of residuals
  plotCorRes(cor_mat = r, pop = as.vector(combined_ordered_breeding[, 1]), ord = ord, 
             title = paste("Evaluation of 1000G admixture proportions with K =", K), 
             max_z = 0.1, min_z = -0.1, rotatelabpop = 90, 
             cex.lab = 0.5, cex.legend = 1)
  
  # Close the PDF file
  dev.off()
}

# Apply the function for K values from 1 to 20
for (K in 1:20) {
  read_and_plot_evalAdmix(K)
}
