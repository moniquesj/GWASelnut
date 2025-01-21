library(vcfR)
library(factoextra)
library(stringr)
library(ggplot2)
library(ggsci)
library(readr)
library(viridis)
library(adegenet)
library(patchwork)

#set wd
setwd("D:/OneDrive - Scuola Superiore Sant'Anna/BioLabs/Projects/Hazelnut/Data analysis/Second_sequencing/datasets")

# Data preprocessing
set.seed(123)
myvcf <- read.vcfR ("population.snps.hz.tombul.filtered.pruned.leaf.vcf")
genind <- vcfR2genind(myvcf)
genind #141 individuals; 16,378 loci; 32,756 alleles 
summary(genind@loc.n.all) # all loci are biallelic
genind_scale <- scaleGen(genind, NA.method="mean")

# PCA
geno.pca <-dudi.pca(genind_scale,cent=TRUE,scale=TRUE,scannf=FALSE, nf=3)
geno.pca_scores <- data.frame(geno.pca$li)# order scores
geno.pca_scores <- data.frame(scale(geno.pca_scores)) #scale scores


# DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS
grp <- find.clusters(genind)
#I chose all PCS and 3 clusters.

# Running the DAPC with 3 clusters
dapc1 <- dapc(genind, grp$grp, n.pca = 200, n.da = 3)

## Tidying the dataframe
assig <- as.data.frame(grp$grp)#in which cluster is each sample?
geno <- rownames(assig)#retrieves the row names of assig and assigns them to geno
rownames(assig) <- NULL #row names are empty of assign
assig <- cbind(geno,assig) #combine the geno vector (containing the row names) and the assig by column binding them together
colnames(assig)<-c("geno", "Cluster")

geno1 <- rownames(geno.pca_scores) #position of each sample in the plot according to the 3 clusters
rownames(geno.pca_scores) <- NULL
geno.pca_scores<- cbind(geno1,geno.pca_scores)
geno.pca_scores_assig <- merge(x= assig, y= geno.pca_scores, by.x = "geno", by.y= "geno1" )

#saving image of the dapc
save.image(file="dapc.Rdata")

#####################################################
##Plotting PCA by varieties and provenance ##########
#####################################################
load(file="dapc.Rdata")

#loading colours again
#colors by country
myCol <- c(
  "Turkey" = "darkgreen",
  "Spain" = "darkturquoise",
  "Italy" = "red",
  "France" = "yellowgreen",
  "United States" = "blue",
  "Greece" = "darkgray",
  "England" = "magenta",
  "Romania" = "darkred",
  "Breeding" = "darksalmon",
  "Unknown" = "lightpink",
  "Germany" = "purple",
  "Former Yugoslavia" = "bisque4",
  "Belgium" = "black",
  "United Kingdom" = "lightblue",
  "Azerbaijan" = "darkkhaki",
  "Georgia" = "yellow",
  "Albania" = "orange",
  "Wild type" = "green"
)

#reading the annotations file
pheno <- read_delim("Simplified_passport_hazelnut_13012025.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)
pheno$Probable.provenance.by.Ferrero <- factor(pheno$Probable.provenance.by.Ferrero) #make Probable provenance as a factor
sample_names <- colnames(myvcf@gt)[-1]
pheno <- pheno[pheno$All.VCF.numbers %in% sample_names, ]

#preparing the files
varieties <- pheno[, c("All.VCF.numbers", "Variety.no.accent", "Probable.provenance.by.Ferrero", "Orchard.country")]
varieties$Variety.no.accent <- as.factor(varieties$Variety.no.accent)
colnames(varieties)[1] <- "geno1" #VCF.number was substituted by geno1
colnames(varieties)[2] <- "Variety"
colnames(varieties)[3] <- "Provenance"
colnames(varieties)[4] <- "Orchard"
varieties_dapc <- varieties

#merging passport information to the cluster and axis information
biostatus <- merge(geno.pca_scores, varieties, by="geno1", all.x=TRUE) 
#biostatus contains the pca scores, the variety name, and the provenances
dapcstatus <- merge(geno.pca_scores_assig, varieties_dapc, by.x="geno", by.y="geno1", all.x=TRUE)
#adding the info of clusters

#write.table(dapcstatus, file="DAPC.clusters.txt", quote=F, row.names=F, sep = "\t")

#plotting the DAPC 1 and 2
dapc1and2 <- ggplot(dapcstatus, aes(x= Axis1, y = Axis2, color=Provenance, shape = Cluster)) +
  scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  #scale_color_nejm() +
  geom_point(size=2) +
  coord_cartesian()+
  labs(x= "PC1 (6.09%)",
       y=  "PC2 (4.38%)", 
       size = "") +
  theme_light() +
  stat_ellipse(data = dapcstatus, aes(x = Axis1, y = Axis2, color = factor(Cluster)),
                                 type = "norm", level = 0.95, linetype = 2)  

dapc1and2


#plotting the DAPC 1 and 3
dapc1and3 <- ggplot(dapcstatus, aes(x= Axis1, y = Axis3, color=Provenance, shape = Cluster)) +
  scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  #scale_color_nejm() +
  geom_point(size=2) +
  coord_cartesian()+
  labs(x= "PC1 (6.09%)",
       y=  "PC3 (3.39%)", 
       size = "") +
  theme_light() + 
  stat_ellipse(data = dapcstatus, aes(x = Axis1, y = Axis3, color = factor(Cluster)),
               type = "norm", level = 0.95, linetype = 2)

dapc1and3

#saving all the images
save_path <- "/output/dapc"
ggsave(pcavar1and2, file = file.path(save_path, "PCA1and2.provenance.jpg"), height = 8, width = 8)
ggsave(pcavar1and3, file = file.path(save_path, "PCA1and3.provenance.jpg"), height = 8, width = 8)
ggsave(dapc1and2, file = file.path(save_path, "DAPC1and2.provenance.jpg"), height = 8, width = 8)
ggsave(dapc1and3, file = file.path(save_path, "DAPC1and3.provenance.jpg"), height = 8, width = 8)


#making a figure with dapcs 
dapc1and2.no.guides <- dapc1and2 + guides(color = 'none', shape = 'none')
dapcs.all <- (dapc1and2.no.guides | dapc1and3)
ggsave(dapcs.all, file = file.path(save_path, "figure.dapc.all.600dpi_wtseparated.jpg"), height=7, width =10, dpi=600)

#making a figure with a big legend
dapc1and2.no.guides <- dapc1and2 + guides(color = 'none', shape = 'none')
dapcs.all <- (dapc1and2.no.guides | dapc1and3)
ggsave(dapcs.all, file = file.path(save_path, "figure.dapc.all.600dpi.small_wtseparated.jpg"), height=4, width =8, dpi=600)

#saving it in svg
#making a figure with a big legend
dapc1and2.no.guides <- dapc1and2 + guides(color = 'none', shape = 'none')
dapcs.all <- (dapc1and2.no.guides | dapc1and3)
ggsave(dapcs.all, file = file.path(save_path, "figure.dapc.all.600dpi.small_wtseparated.svg"), height=4, width =8, dpi=600)

#making a small figure with big legend
dapc1and2 <- ggplot(dapcstatus, aes(x= Axis1, y = Axis2, color=Provenance, shape = Cluster)) +
  scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  #scale_color_nejm() +
  geom_point(size=4) +
  coord_cartesian()+
  labs(x= "PC1 (6.09%)",
       y=  "PC2 (4.38%)", 
       size = "") +
  theme_light() +
  theme(legend.text = element_text(size = 15),  # Adjust the text size
        legend.title = element_text(size = 17),  # Adjust the title size
        legend.key.size = unit(1.5, 'lines')) +  # Adjust the key size
  guides(shape = guide_legend(override.aes = list(size = 5))) +  # Adjust the shape size in legend 
  stat_ellipse(data = dapcstatus, aes(x = Axis1, y = Axis2, color = factor(Cluster)),
               type = "norm", level = 0.95, linetype = 2)  

dapc1and3 <- ggplot(dapcstatus, aes(x= Axis1, y = Axis3, color=Provenance, shape = Cluster)) +
  scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  #scale_color_nejm() +
  geom_point(size=4) +
  coord_cartesian()+
  labs(x= "PC1 (6.09%)",
       y=  "PC3 (3.39%)", 
       size = "") +
  theme_light() + 
  theme(legend.text = element_text(size = 15),  # Adjust the text size
        legend.title = element_text(size = 17),  # Adjust the title size
        legend.key.size = unit(1.5, 'lines')) +  # Adjust the key size
  guides(shape = guide_legend(override.aes = list(size = 5))) +  # Adjust the shape size in legend
  stat_ellipse(data = dapcstatus, aes(x = Axis1, y = Axis3, color = factor(Cluster)),
               type = "norm", level = 0.95, linetype = 2)


dapc1and2.no.guides <- dapc1and2 + guides(color = 'none', shape = 'none')
dapcs.all <- (dapc1and2.no.guides | dapc1and3)
ggsave(dapcs.all, file = file.path(save_path, "big.legend.dapc.600dpi_wtseparated.jpg"), height=10, width =8, dpi=600)

#saving it in svg
ggsave(dapcs.all, file = file.path(save_path, "big.legend.dapc.600dpi_wtseparated.svg"), height=10, width =8, dpi=600)

# Plot of the BIC values
bic_values <- grp$Kstat
bic_values <- as.data.frame(bic_values)
K <- rownames(bic_values)
rownames(bic_values) <- NULL
bic_values <- cbind(K,bic_values)
colnames(bic_values)[2] <- "bic"
bic_values$K <- sub("K=", "", bic_values$K) #removing the "K=" prefix
bic_values$K <- as.numeric(bic_values$K)

# Plotting
bic.plot <- ggplot(data=bic_values, mapping = aes(x=K, y=bic)) +
  geom_line()+
  ylab("BIC")+
  xlab("K")+
  scale_x_continuous(breaks = seq(min(bic_values$K), max(bic_values$K), by = 2)) +
  ggtitle("Value of Bayesian information criterion (BIC) \n versus number of clusters (K)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 18, hjust = 0.5),        # Increase title size
        axis.title.x = element_text(size = 16),                   # Increase x-axis label size
        axis.title.y = element_text(size = 16),                   # Increase y-axis label size
        axis.text = element_text(size = 14)                       # Increase tick label size
  )

bic.plot

# Save the plot to a specified directory
ggsave(filename = file.path(save_path, "bic-values-plot-dapc.png"), 
       plot = bic.plot, dpi = 600, height=5, width=6)

ggsave(filename = file.path(save_path, "bic-values-plot-dapc.svg"), 
       plot = bic.plot, dpi = 600, height=5, width=6)


