# Admixture maps with pie charts

# The goal of this script is to plot the results from ADMIXTURE with K=5 (because it is the one with lowest CV error) in a map
# Load necessary libraries
library(tidyr)
library(tidyverse)
library(ggplot2)
library(mapmixture)
library(gridExtra)
library(plotrix)
library(maps)
library(readr)
library(sf)

# Set wd
setwd("output/admixture")

# Import coordinates and the dataset from admixture
coords = read.csv("coordinates.grouped.without.breeding.wt.unknown.two.european.clusters.csv")

combined_ordered_breeding <- read_delim("Combined_ordered_breeding_K5.txt", 
                              delim = "\t", escape_double = FALSE, 
                              col_types = cols(Cluster1 = col_number(), 
                                               Cluster2 = col_number(), Cluster3 = col_number(), 
                                               Cluster4 = col_number(), Cluster5 = col_number()), 
                                              trim_ws = TRUE)
# Adjusting names of columns
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "sample"] <- "Ind"
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "provenance"] <- "Site"
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "V1"] <- "Cluster1"
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "V2"] <- "Cluster2"
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "V3"] <- "Cluster3"
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "V4"] <- "Cluster4"
names(combined_ordered_breeding)[names(combined_ordered_breeding) == "V5"] <- "Cluster5"

# Make Site as the first column 
combined_ordered_breeding <- combined_ordered_breeding %>%
  select(Site, Ind, everything())

# Drop WT, Breeding, and Unknown, as we do not have their provenance information
dataset_map <- subset(combined_ordered_breeding, !(Site %in% c("Wild type", "Unknown", "Breeding")))

# Adjusting the clusters
# Eurasia - Azerbaijan (3), Georgia (11), Turkey (16)
# Northern Europe - France (5), Germany (1), England (3)
# Central Europe - Greece (1), Romania (1), Former Jugoslavia (1)
# Italy
# Spain
# USA

dataset_map_grouped <- dataset_map %>%
  mutate(Site = ifelse(Site %in% c("Azerbaijan", "Georgia", "Turkey"), "Eurasia", Site))

dataset_map_grouped <- dataset_map_grouped %>%
  mutate(Site = ifelse(Site %in% c("France", "Germany", "England"), "Northern Europe", Site))

dataset_map_grouped <- dataset_map_grouped %>%
  mutate(Site = ifelse(Site %in% c("Greece", "Romania", "Former Jugoslavia"), "Central Europe", Site))

# Order both files alphabetically by site
coords = coords[order(coords$Code), ] 
dataset_map_grouped = dataset_map_grouped[order(dataset_map_grouped$Site), ] 

# Plotting the map with bigger pie size and legend
map1 <- mapmixture(dataset_map_grouped, 
                   coords, 
                   crs=4326, # My coordinates use the WGS84 system, therefore crs=4326
                   cluster_cols = c("#FF0000", "#CCFF00", "#00FF66", "#0066FF", "#CC00FF"),
                   cluster_names = c("Ancestral Cluster 1","Ancestral Cluster 2", "Ancestral Cluster 3", "Ancestral Cluster 4", "Ancestral Cluster 5"),
                   pie_size = 5,
                   axis_title_size = 20,
                   axis_text_size = 18) +
  # Adjust theme options
  theme(
    legend.position = "top",
    plot.margin = margin(l = 5, r = 5),
    legend.text = element_text(size = 20)
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 18, alpha = 1))) 

ggsave("admixture.map.grouped.two.european.clusters_bigger.svg",map1, width=15, height=10, dpi=600)

# The labels will be added manually

# Saving again just to make a big legend
map1 <- mapmixture(dataset_map_grouped, 
                   coords, 
                   crs=4326, # My coordinates use the WGS84 system, therefore crs=4326
                   cluster_cols = c("#FF0000", "#CCFF00", "#00FF66", "#0066FF", "#CC00FF"),
                   cluster_names = c("Ancestral Cluster 1","Ancestral Cluster 2", "Ancestral Cluster 3", "Ancestral Cluster 4", "Ancestral Cluster 5"),
                   pie_size = 5,
                   axis_title_size = 20,
                   axis_text_size = 18) +
  # Adjust theme options
  theme(
    legend.position = "top",
    plot.margin = margin(l = 5, r = 5),
    legend.text = element_text(size = 20)
  )+
  # Adjust the size of the legend keys
  guides(fill = guide_legend(override.aes = list(size = 18, alpha = 1))) 

ggsave("admixture.map.grouped.two.european.clusters_bigger_legend_ancestral_cluster.svg",map1, width=35, height=10, dpi=600)
