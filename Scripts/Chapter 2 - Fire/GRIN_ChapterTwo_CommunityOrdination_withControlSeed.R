################################################################################
################################################################################
#########################     GRIN - Fire         ##############################
#########################    NMDS - Community     ##############################
#########################    with Tw Control      ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2021 - 2023         ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "labdsv", 
                      "pairwise.adonis", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(devtools)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

##########################     Read in 2021-2023 Data       ####################

GRIN = read.csv("Data/GRIN - 2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)

str(GRIN)
summary(GRIN)

# Remove Tilling Treatment # 
GRIN = filter(GRIN, Treatment != 'Tsp')
GRIN = filter(GRIN, Treatment != 'C')

# Separate Years #
GRIN_21 = filter(GRIN, Year == 1)
GRIN_22 = filter(GRIN, Year == 2)
GRIN_23 = filter(GRIN, Year == 3)

########################### 2023 Data ##########################################

# Create species pivot table #
Spp = dplyr::select(GRIN_23, ID, Species, Coverage) %>% matrify()

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_23, ID, Treatment, Fire, T.F) %>% summarise()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot three dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(Spp, distance = "bray", k=3)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN_23, Species, Group) %>% summarise()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- 
  as.data.frame(vegan:::scores.metaMDS(MDS, display = c("species")))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #                                                       )  
species.scores$Group <- species_groups$Group

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment, 
                  Fire = Treat$Fire, T.F = Treat$T.F, Plot = Treat$ID)

# NMDS Graphs
NMDS_plot = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = T.F),
             alpha = 0.7, size = 3, shape = 21) +
  # geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species)) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, fill = T.F, color = T.F), 
               linetype = "solid", show.legend = T, alpha = 0.15) +
  #scale_fill_manual(labels=c('Control', 'Late-Spring', 'Winter'),
  #                  values=c("#FF3399", "#117733", "#3366FF")) +
  #scale_color_manual(labels=c('Control', 'Late-Spring', 'Winter'),
  #                  values=c("#FF3399", "#117733", "#3366FF")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        legend.text=element_text(face = "bold", color = "black")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_plot

ggsave("Figures/Chapter 2 - Fire/2023_NMDS_Tw.png", 
       width = 10, height = 7)

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Spp ~ NMDS$T.F, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp ~ T.F, data = NMDS)
pairwise.adonis


##########################     2022 Data       #################################

# Create species pivot table #
Spp = dplyr::select(GRIN_22, ID, Species, Coverage) %>% matrify()

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_22, ID, Treatment, Fire, T.F) %>% summarise()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot two dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(Spp, distance = "bray", trymax = 500, maxit = 999, k=4, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN_22, Species, Group) %>% summarise()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- 
  as.data.frame(vegan:::scores.metaMDS(MDS, display = c("species")))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #
species.scores$Group <- species_groups$Group

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Fire = Treat$Fire, T.F = Treat$T.F, 
                  Treatment = Treat$Treatment, Plot = Treat$ID)

# NMDS Graphs
NMDS_plot = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = T.F),
             alpha = 0.7, size = 3, shape = 21) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, fill = T.F, color = T.F), 
               linetype = "solid", show.legend = T, alpha = 0.15) +
  #scale_fill_manual(labels=c('Control', 'Late-Spring', 'Winter'),
  #                  values=c("#FF3399", "#117733", "#3366FF")) +
  #scale_color_manual(labels=c('Control', 'Late-Spring', 'Winter'),
  #                   values=c("#FF3399", "#117733", "#3366FF")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        legend.text=element_text(face = "bold", color = "black")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_plot

ggsave("Figures/Chapter 2 - Fire/2022_NMDS_Tw.png", 
       width = 10, height = 7)

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Spp ~ NMDS$T.F, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp ~ T.F, data = NMDS)
pairwise.adonis

##########################     2021 Data       #################################

# Create species pivot table #
Spp = dplyr::select(GRIN_21, ID, Species, Coverage) %>% matrify()

#Create grouped treatment/environment table and summaries to fit species table #
Treat = group_by(GRIN_21, ID, treatment, Fire, T.F) %>% summarise()

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),
           replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot two dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(Spp, distance = "bray", trymax = 500, maxit = 999, k=3, 
              trace = F, autotransform = FALSE, wascores = TRUE)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Create a frame for functional groups alongside species for NMDS graph #
species_groups = group_by(GRIN, Species, Group) %>% summarise()

# Extract  species scores & convert to a data.frame for NMDS graph #
species.scores <- 
  as.data.frame(vegan:::scores.metaMDS(MDS, display = c("species")))

# create a column of species, from the row names of species.scores  #                                                            )  
species.scores$species <- rownames(species.scores)

# create a column for functional groups for NMDS graph #                                                       )  
species.scores$Group <- species_groups$Group

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment, 
                  Fire = Treat$Fire, T.F = Treat$T.F, Plot = Treat$ID)

# NMDS Graphs
NMDS_plot = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, fill = T.F),
             alpha = 0.7, size = 3, shape = 21) +
  # geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species)) +
  stat_ellipse(geom = "polygon", data = NMDS, 
               aes(x = MDS.MDS1, y = MDS.MDS2, fill = T.F, color = T.F), 
               linetype = "solid", show.legend = T, alpha = 0.15) +
  #scale_fill_manual(labels=c('Control', 'Late-Spring', 'Winter'),
  #                  values=c("#FF3399", "#117733", "#3366FF")) +
  #scale_color_manual(labels=c('Control', 'Late-Spring', 'Winter'),
  #                   values=c("#FF3399", "#117733", "#3366FF")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        legend.text=element_text(face = "bold", color = "black")) +
  labs(x = "MDS1", y = "MDS2", color = "Fire Treatment", 
       fill = "Fire Treatment")
NMDS_plot

ggsave("Figures/Chapter 2 - Fire/2021_NMDS_Tw.png", 
       width = 10, height = 7)

# Perform adonis to test the significance of treatments#
adon.results <- adonis2(Spp ~ NMDS$T.F, method="bray",perm=999)
print(adon.results)
pairwise.adonis<-pairwise.adonis2(Spp ~ T.F, data = NMDS)
pairwise.adonis

