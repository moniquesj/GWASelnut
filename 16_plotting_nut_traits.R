#In this script, I will take a look at the seed traits data
library(readr)
library(dplyr)
library(ggplot2)
library(readr)
library(ggsci)
library(ggpubr)
library(car)
library(rstatix)
library(PerformanceAnalytics)
library(Hmisc)
library(gridExtra)
library(GAPIT)
library(tidyr)
library(ggpubr)
library(bestNormalize)
library(purrr)

#importing the list of samples effectively used for the GWAS of the nut morphology
setwd("/gwas_all_traits_clean_output/samples for each gwas")
samples_nut_list <- read_delim("samples_nut_traits_gwas.txt", 
                                    delim = "\t", escape_double = FALSE, 
                                    col_names = FALSE, trim_ws = TRUE)
colnames(samples_nut_list)[1] <- "All.VCF.numbers"


#importing the passport file
setwd("/Second_sequencing/datasets")
passp <- read_delim("Simplified_passport_hazelnut_13012025.txt",
                   delim = "\t", escape_double = FALSE,
                   trim_ws = TRUE)

#considering only leaves samples for GWAS
passp.leaf <- subset(passp, Leaf == 1)

#keeping only the meaningful columns
columns_to_keep <- c(1, 7, 24, 26, 90:93)
filtered.data <- passp.leaf[, columns_to_keep]

#removing all samples that were not used for the GWAS
filtered.data <- filtered.data %>%
  filter(All.VCF.numbers %in% samples_nut_list$All.VCF.numbers)

# Are the traits correlated?
corr.data <- filtered.data[, 5:8]
corr.data <- as.data.frame(corr.data)
corr.plot <- chart.Correlation(corr.data[,1:4],pch=19, method = 'pearson')

# Save as SVG
svg(file = file.path(save_path, "corr_plot_only_nut_traits_with_associated_snps.svg"), height = 8, width = 10)
chart.Correlation(corr.data[,1:4], pch = 19, method = 'pearson')
dev.off()

# Save as PNG
png(file = file.path(save_path, "corr_plot_only_nut_traits_with_associated_snps.png"), height = 8, width = 10, units = "in", res = 400)
chart.Correlation(corr.data[,1:4], pch = 19, method = 'pearson')
dev.off()


# Grouping the samples according to the groups used for the ADMIXTURE map
filtered.data <- filtered.data %>%
  mutate(Provenance.groups = case_when(
    Probable.provenance.by.Ferrero %in% c("Azerbaijan", "Georgia", "Turkey") ~ "Eurasia",
    Probable.provenance.by.Ferrero %in% c("France", "Germany", "England") ~ "Central Europe",
    Probable.provenance.by.Ferrero %in% c("Greece", "Romania", "Former Jugoslavia", "Albania") ~ "Eastern Europe",
    TRUE ~ Probable.provenance.by.Ferrero  # Retains original value if not matching
  ))

# I'm adding the only sample from Albania into the Eastern Europe group


#################################################
######## SHAPIRO TEST FOR ALL VARIABLES #########
#################################################
# Nut size
filtered.data %>%
  group_by(Provenance.groups) %>%
  filter(n() >= 3) %>%  # keep groups with at least 3 samples, otherwise Shapiro test cannot be performed
  summarise(statistic = shapiro.test(Nut.average.size)$statistic)

# Nut perimeter
filtered.data %>%
  group_by(Provenance.groups) %>%
  filter(n() >= 3) %>%  # keep groups with at least 3 samples, otherwise Shapiro test cannot be performed
  summarise(statistic = shapiro.test(Nut.perimeter)$statistic)

# Nut maximum caliber
filtered.data %>%
  group_by(Provenance.groups) %>%
  filter(n() >= 3) %>%  # keep groups with at least 3 samples, otherwise Shapiro test cannot be performed
  summarise(statistic = shapiro.test(Nut.maximum.caliber)$statistic)

# Nut minimum caliber
filtered.data %>%
  group_by(Provenance.groups) %>%
  filter(n() >= 3) %>%  # keep groups with at least 3 samples, otherwise Shapiro test cannot be performed
  summarise(statistic = shapiro.test(Nut.minimum.caliber)$statistic)

# All traits are normal within each provenance that has at least 3 samples

#################################################
######## LEVENE TEST FOR ALL VARIABLES #########
#################################################
#testing homoscedasticity (are the variances the same within each provenance?)
filtered.data$Provenance.groups <- as.factor(filtered.data$Provenance.groups)
leveneTest(Nut.average.size ~ Provenance.groups, data = filtered.data)
leveneTest(Nut.perimeter ~ Provenance.groups, data = filtered.data)
leveneTest(Nut.maximum.caliber ~ Provenance.groups, data = filtered.data)
leveneTest(Nut.minimum.caliber ~ Provenance.groups, data = filtered.data)

# All variances are the same within each provenance, except for the Nut.maximum.caliber and Nut.average.size

#################################################
######## ANOVAS #################################
#################################################
# Retaining only samples from  provenances that have at least 3 samples belonging to them
filtered.data.min.provenances <- filtered.data %>%
  group_by(Provenance.groups) %>%
  filter(n() >= 3) %>%
  ungroup()

# Performing the ANOVAs
anova.nut.perimeter <- aov(Nut.perimeter ~ Provenance.groups, data = filtered.data.min.provenances) %>%
  tukey_hsd() 
anova.nut.perimeter 

# Performing the Kruskal-Walis test for the traits that are not homoschedastisc
kruskal.nut.average.size <- kruskal.test(Nut.average.size ~ Provenance.groups, data = filtered.data.min.provenances)
kruskal.nut.average.size # No differences between groups, therefore no need to perform Dunn's test

kruskal.nut.max.caliber <- kruskal.test(Nut.maximum.caliber ~ Provenance.groups, data = filtered.data.min.provenances)
kruskal.nut.max.caliber # No differences between groups, therefore no need to perform Dunn's test

kruskal.nut.min.caliber <- kruskal.test(Nut.minimum.caliber ~ Provenance.groups, data = filtered.data.min.provenances)
kruskal.nut.min.caliber # No differences between groups, therefore no need to perform Dunn's test

# No differences between the groups in any of the traits

#################################################
######## BOXPLOTS ###############################
#################################################
# setting the colours
myCol <- c(
  "Spain" = "darkturquoise",
  "Italy" = "red",
  "United States" = "blue",
  "Breeding" = "darksalmon",
  "Unknown" = "lightpink",
  "Eurasia" = "goldenrod",
  "Eastern Europe" ="deepskyblue",
  "Central Europe" = "mediumseagreen"
)


# setting the y breaks
nut.size.y.breaks <- seq(min(filtered.data$Nut.average.size), max(filtered.data$Nut.average.size), by = 32) 
nut.perimeter.y.breaks <- seq(min(filtered.data$Nut.perimeter), max(filtered.data$Nut.perimeter), by = 8) 
nut.max.caliber.y.breaks <- seq(min(filtered.data$Nut.maximum.caliber), max(filtered.data$Nut.maximum.caliber), by = 3) 
nut.min.caliber.y.breaks <- seq(min(filtered.data$Nut.minimum.caliber), max(filtered.data$Nut.minimum.caliber), by = 2) 

# Boxplots
# Nut size
size <- ggboxplot(filtered.data, x = "Provenance.groups", y = "Nut.average.size", fill="Provenance.groups") +
  ylab (bquote(bold("Nut size (mm²)"))) +
  ggtitle(bquote(bold("Nut size across provenances"))) +
  xlab (bquote(bold("Provenance"))) +
  border("black") +
  geom_point(position = "jitter") +
  theme(legend.position = "none")+
  stat_compare_means(method = "kruskal.test", label.y = 225)+
  scale_fill_manual(values = myCol)+
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = nut.size.y.breaks, labels = scales::number_format(scale = 1, accuracy = 1))

size

# Nut perimeter
perimeter <- ggboxplot(filtered.data, x = "Provenance.groups", y = "Nut.perimeter", fill="Provenance.groups") +
  ylab (bquote(bold("Nut perimeter (mm)"))) +
  ggtitle(bquote(bold("Nut perimeter across provenances"))) +
  xlab (bquote(bold("Provenance"))) +
  border("black") +
  geom_point(position = "jitter") +
  theme(legend.position = "none")+
  stat_compare_means(method = "anova", label.y = 70)+
  scale_fill_manual(values = myCol)+
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = nut.perimeter.y.breaks, labels = scales::number_format(scale = 1, accuracy = 1))

perimeter

# Nut maximum caliber
max.caliber <- ggboxplot(filtered.data, x = "Provenance.groups", y = "Nut.maximum.caliber", fill="Provenance.groups") +
  ylab (bquote(bold("Nut maximum caliber (mm)"))) +
  ggtitle(bquote(bold("Nut maximum caliber across provenances"))) +
  xlab (bquote(bold("Provenance"))) +
  border("black") +
  geom_point(position = "jitter") +
  theme(legend.position = "none")+
  stat_compare_means(method = "kruskal.test", label.y = 25)+
  scale_fill_manual(values = myCol)+
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = nut.max.caliber.y.breaks, labels = scales::number_format(scale = 1, accuracy = 1))

# Nut minimum caliber
min.caliber <- ggboxplot(filtered.data, x = "Provenance.groups", y = "Nut.minimum.caliber", fill="Provenance.groups") +
  ylab (bquote(bold("Nut minimum caliber (mm)"))) +
  ggtitle(bquote(bold("Nut minimum caliber across provenances"))) +
  xlab (bquote(bold("Provenance"))) +
  border("black") +
  geom_point(position = "jitter") +
  theme(legend.position = "none")+
  stat_compare_means(method = "kruskal.test", label.y = 16)+
  scale_fill_manual(values = myCol)+
  guides(fill=FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = nut.min.caliber.y.breaks, labels = scales::number_format(scale = 1, accuracy = 1))


#creating a figure
figure <- grid.arrange(
  arrangeGrob(size, perimeter,max.caliber, min.caliber, ncol = 2, nrow = 2),  
  heights = c(2),     #Adjust row heights
  widths = c(2))    # Adjust column widths

#saving the figure
setwd("/output/gwas_nut_traits")
ggsave(figure, file="boxplots_nut_traits_admixture_groups_dpi600.svg", height=10, width=11.5, dpi = 600)
ggsave(figure, file="boxplots_nut_traits_admixture_groups_dpi400.png", height=10, width=11.5, dpi = 400)

# getting means and sd
mean(filtered.data$Nut.average.size)
sd(filtered.data$Nut.average.size)
mean(filtered.data$Nut.perimeter)
sd(filtered.data$Nut.perimeter)
mean(filtered.data$Nut.maximum.caliber)
sd(filtered.data$Nut.maximum.caliber)
mean(filtered.data$Nut.minimum.caliber)
sd(filtered.data$Nut.minimum.caliber)
