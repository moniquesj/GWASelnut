# Load necessary libraries
library(tidyr)
library(tidyverse)
library(ggplot2)

# Set wd
setwd("/output/admixture")

# Plotting the CV error
cv <- read_delim("cv-error-admixture.txt", 
                 delim = "\t", escape_double = FALSE, 
                 col_types = cols(K = col_number(), `CV_error` = col_number()), 
                 trim_ws = TRUE)

cv.plot <- ggplot(data=cv, mapping = aes(x=K, y=CV_error)) +
  geom_line()+
  ylab("CV error")+
  xlab("K")+
  scale_x_continuous(breaks = seq(min(cv$K), max(cv$K), by = 2))+
  ggtitle("Value of cross-validation (CV) error \n versus number of clusters (K)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))

cv.plot

ggsave("cv-error-plot-admixture.png", plot = cv.plot, dpi = 600, height=5, width=5)
ggsave("cv-error-plot-admixture.svg", plot = cv.plot, dpi = 600, height=5, width=5)

# Read the q-estimates file 
q_file <- read.table("population.snps.hz.tombul.filtered.pruned.leaf.for.plink.admixture.5.Q", header = FALSE)

# Read the meta-file with sample names and groups 
meta_file <- read_delim("Simplified_passport_hazelnut_for_plink_admixture.txt", delim = "\t", 
                    escape_double = FALSE, trim_ws = TRUE)
meta_file <- as.data.frame(meta_file)

# Adjusting names of columns of meta file
names(meta_file)[names(meta_file) == "All.VCF.numbers.ADMIXTURE"] <- "sample"
names(meta_file)[names(meta_file) == "Population.ADMIXTURE"] <- "provenance"

# Read the ped file with sample names in the correct order 
ped_file <- read.table("population.snps.hz.tombul.filtered.pruned.leaf.for.plink.admixture.ped", header = FALSE)

# Filter the meta-file to keep only samples that are in the ped file
filtered_meta <- meta_file[meta_file$sample %in% ped_file$V1, ]

# Order the meta-file according to the sample names in the ped file
ordered_meta <- filtered_meta[match(ped_file$V1, filtered_meta$sample), ]

# Combine meta data and q data
combined <- cbind(ordered_meta, q_file)

# Order according to provenance categories
combined_ordered <- combined[order(combined$provenance), ]

# Reorder the dataset, placing rows with provenance "Breeding" at the end
combined_ordered_breeding <- rbind(
  combined_ordered[combined_ordered$provenance != "Breeding", ],
  combined_ordered[combined_ordered$provenance == "Breeding", ]
)

# Determine the mean row number that is in the middle of each provenance category
xlabels <- aggregate(1:nrow(combined_ordered_breeding),
                     by = list(combined_ordered_breeding[, "provenance"]),
                     FUN = mean)
xlabels


# Determine the highest row for each provenance category
sampleEdges <- aggregate(1:nrow(combined_ordered_breeding),
                         by = list(combined_ordered_breeding[, "provenance"]), 
                         FUN = max)
sampleEdges

# Create the bar plot for the figure of the manuscript
# Axis names and the Greece label will be added manually on Inkscape

svg("q-estimates-K5-b2000_big_label.svg", width = 35, height = 10)
par(mar = c(25, 10, 4, 2) + 0.1)  # Increase bottom and left margins
par(cex.lab = 3, cex.axis = 3, cex.main = 3)  # Increase label, axis, and title sizes

# Set las = 2 for vertical tick labels on the x-axis
par(las = 2)  # This rotates only the tick labels for the x-axis

# Generate the barplot
barplot(t(as.matrix(combined_ordered_breeding[, -1:-2])), col = rainbow(5), 
        space = 0, ylab = "", xlab = "",  # We'll manually add axis labels with mtext
        border = NA, axisnames = FALSE)

abline(v = sampleEdges$x, lwd = 2) # Add vertical lines
axis(1, at = xlabels$x - 0.5, labels = xlabels$Group.1) # Add x-axis tick labels (Group names, vertical due to las = 2)
title(main = "K=5", col.main = "black")
dev.off()



