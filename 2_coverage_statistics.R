#load libraries
library(tidyverse)
library(readr)
library(ggplot2)

# I followed the section 'Variant mean depth' from this tutorial
# https://speciationgenomics.github.io/filtering_vcfs/

# Set wd
setwd("/analysis/output-vcftools")

# Checking variant depth
var_depth <- read_delim("output-vcftools.ldepth.mean", delim = "\t",
                        col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

plot.var_depth <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
plot.var_depth + theme_light()

#This plot is a bit misleading because clearly, there are very few variants 
#with extremely high coverage indeed. 
#Let’s take a closer at the mean depth
summary(var_depth$mean_depth)
#Median is 23 and mean is 34

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 3.103   10.938   23.083   33.997   50.927 1713.900 

#Redraw the plot with xlim
plot.var_depth + theme_light() + xlim(0, 100)

# This gives a better idea of the distribution.
# We could set our minimum coverage at the 5 and 95% quantiles
# but we should keep in mind that the more reads that cover a site,
# the higher confidence our basecall is. 10x is a good rule of thumb as a
# minimum cutoff for read depth, although if we wanted to be conservative,
# we could go with 15x.

# Minimum cutoff: 15x is a good rule of thumb as a minimum cutoff for read depth
# Maximum cutoff: Usually a good rule of thumb is something the mean depth x 2 - so in 
# this case we could set our maximum depth at 68x, since the mean is 33.99.

# Therefore, I will set the minimum depth to 15X and the maximum to 68X