library(vcfR)
library(readr)
library(adegenet)
library(ape)
library(ggplot2)
library(ggtree)
library(randomcoloR)
library(dplyr)
library(poppr)

set.seed(123) 

#setting wd
setwd("Second_sequencing/datasets")

#opening the vcf file
vcf <- read.vcfR("population.snps.hz.tombul.filtered.pruned.leaf.vcf")

#load the annotations file
passp <- read_delim("Simplified_passport_hazelnut_13012025.txt", delim = "\t",
                            escape_double = FALSE, trim_ws = TRUE)
passp <- as.data.frame(passp)

# Read the sample names from the VCF file
sample_names <- colnames(vcf@gt)[-1]

#filter the passp to keep only the samples present in the VCF file
filtered_passp <- passp[passp$All.VCF.numbers %in% sample_names, ]

#convert vcf to genlight object
hz.gl <- vcfR2genlight(vcf)

#In order to specify the population, I am adding the Probable provenance column from the
#pop.data data frame to the pop slot of the genlight object
pop(hz.gl) <- filtered_passp$Probable.provenance.by.Ferrero

# # Checking unique values in the Probable.provenance column
unique_countries <- unique(filtered_passp$Probable.provenance.by.Ferrero)

# Create an empty vector to store the colors
colors <- vector(length = nrow(filtered_passp))

# Assigning colors based on country
for (i in 1:length(unique_countries)) {
  country <- unique_countries[i]

if (country == "Turkey") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "darkgreen"
} else if (country == "Spain") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "darkturquoise"
} else if (country == "Italy") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "red"
} else if (country == "France") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "yellowgreen"
} else if (country == "United States") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "blue"
} else if (country == "Greece") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "darkgray"
} else if (country == "England") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "magenta"
} else if (country == "Romania") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "darkred"
} else if (country == "Breeding") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "darksalmon"
} else if (country == "Unknown") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "lightpink"
} else if (country == "Germany") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "purple"
} else if (country == "Former Yugoslavia") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "bisque4"
} else if (country == "Belgium") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "black"
} else if (country == "United Kingdom") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "lightblue"
} else if (country == "Azerbaijan") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "darkkhaki"
} else if (country == "Georgia") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "yellow"
} else if (country == "Albania") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "orange"
} else if (country == "Wild type") {
  colors[filtered_passp$Probable.provenance.by.Ferrero == country] <- "green"
}
}

# Add the colors column to the data frame
filtered_passp$colors <- colors

#setting the ploidy
ploidy(hz.gl) <- 2

#boostrapt analysis with poppr
tree.boot <- aboot(hz.gl, tree = "nj", distance = "bitwise.dist", sample = 2000, showtree = TRUE)

# #saving my objects in a file 
save.image(file = "tree.R.data")

#opening the file from the unrooted tree
load("tree.R.data")

#choosing the sample to root
sample_to_root <- "261-HZ266-P3b-E09.1" #this is the only Tombul whose provenance and orchard.country are Turkey

#rooting the tree 
rooted_tree <- root(tree.boot, sample_to_root)

#reorder All.VCF.numbers according to the tree
desired_order <- rooted_tree$tip.label
filtered_passp <- filtered_passp %>%
slice(match(desired_order, All.VCF.numbers))

#changing the tip labels to the variety's names
# Retrieve the current tip labels from the tree.boot object
current_labels <- rooted_tree$tip.label
# Get the corresponding Variety names from filtered_passp
new_labels <- filtered_passp$Variety.no.accent
# Replace the tip labels in the tree.bot object
rooted_tree$tip.label <- new_labels[match(current_labels, filtered_passp$All.VCF.numbers)]

# plotting the tree without variety's names
options(ignore.negative.edge=TRUE)
plot.rooted.tree <-  ggtree(rooted_tree, layout="dendrogram")+
  geom_tippoint(color=filtered_passp$colors, size=4.5) 

plot.rooted.tree

#saving tree in png
plot_filename <- file.path(save_path, "tree.leaf.rooted.novarnames.png")
ggsave(filename = plot_filename, plot = plot.rooted.tree, height = 10, width = 15, dpi = 600)

#saving tree in svg
save_path <- "/output/trees"
plot_filename <- file.path(save_path, "tree.leaf.rooted.novarnames.svg")
ggsave(filename = plot_filename, plot = plot.rooted.tree, height = 10, width = 15, dpi = 600)


# Plotting the supplementary tree
#plotting the tree
options(ignore.negative.edge=TRUE)
plot.rooted.tree <-  ggtree(rooted_tree, layout="rectangular")+
  geom_tiplab(size=3.0, color="black", hjust = -0.1)+
  geom_tippoint(color=filtered_passp$colors, size=3.0) +
  geom_nodelab(size=2.0, hjust = 0)

plot.rooted.tree

#saving tree in png 
save_path <- "/output/trees"
plot_filename <- file.path(save_path, "supplem_tree.leaf.rooted.varnames.png")
ggsave(filename = plot_filename, plot = plot.rooted.tree, height = 10, width = 15, dpi = 600)

#saving tree in svg
plot_filename <- file.path(save_path, "supplem_tree.leaf.rooted.varnames.svg")
ggsave(filename = plot_filename, plot = plot.rooted.tree, height = 10, width = 15, dpi = 600)
