rm(list=ls(all=TRUE)) # Clears environment
cat("\014") # Clears history

list.of.packages <- c("ggplot2", "reshape2", "gridExtra", "tidyverse", "grid", "labdsv")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(tidyverse)
library(vegan)
library(grid)
library(labdsv)

set.seed(2)

Cov = read.csv("Seeding - Functional Groups - Coverage.csv")

# Create Grouped Treatment/ Environment Table and Summarise to fit Species Table #

Treat = group_by(Cov, Fire, Tilling, Seeding, Plot)

# Create Species Pivot Table #

spp = subset(Cov, select = -c(Plot, Fire, Tilling, Seeding) )
spp = as.matrix(spp)
spp = as.data.frame(spp)

# Use dissimilarities to create scree plot - attain the number of dimensions for NMDS with least stress #
# Using function that produces a stress vs. dimensionality plot #

# Run Bray-Curtis Dissimilarity #

vegdist = vegdist(spp, method = "bray")
NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", 
       main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(vegdist)

# - Based on scree plot two dimesions will be sufficent for NMDS #

# Run MDS and plot stress using a Shepherd Plot #

MDS = metaMDS(spp, distance = "bray", trymax = 1000, k=2)
MDS$stress
stressplot(MDS) 
goodness(MDS)

# Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Place species and subplot scores from MDS into a dataframe for ggplot #

spp_scores<- data.frame(scores(MDS, "species"))
spp_scores$species <- row.names(spp_scores)
sites_scores <- data.frame(scores(MDS, "sites"))
sites_scores$sites <- row.names(sites_scores)

plot(MDS)
ordiplot(spp_scores)

NMDS = data.frame(MDS = MDS$points, sites_scores = sites_scores, spp_scores = spp_scores, 
                  Plot = Treat$Plot, Fire = Treat$Fire, Tilling = Treat$Tilling, Seeding = Treat$Seeding )

Fire = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, color = Fire, shape = Seeding)) +
  theme_bw() +
  labs(x="MDS1", y="MDS2", title = "NMDS GRIN Plots", color = "Fire", shape = "Seeding") +
  theme(plot.title = element_text(hjust = 0.5)) 
Fire

ggplot() +
  geom_bar(data = Cov, aes(Grass, Woody))
