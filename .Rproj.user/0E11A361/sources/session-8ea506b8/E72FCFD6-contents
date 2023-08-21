#########################   GRIN - Disturbance    ##############################
#########################      Seed Bank          ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2022        ##############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("extrafont", "tidyverse", "vegan", "labdsv")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(extrafont)
library(tidyverse)
library(vegan)
library(labdsv)
#font_import()
#loadfonts(device = "win")

##########################     Read in Data       #########################

GRIN = read.csv("Data/GRIN - SeedBank.csv")
str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Create grouped treatment/environment table and summaries to fit species table #
Treat = group_by(GRIN, ID, Treatment) %>% dplyr::summarize()

GRIN$Plot = NULL
GRIN$Sub_Plot = NULL
GRIN$Treatment = NULL
GRIN = GRIN %>% remove_rownames %>% column_to_rownames(var="ID")

# Use dissimilarities to create scree plot - attain the number of dimensions #
# for NMDS with least stress. Using function that produces a # 
# stress vs. dimensional plot #

NMDS.scree <- function(x) { # x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), 
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", 
       ylab = "Stress", main = "NMDS Stress Plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#NMDS.scree(Spp) 
# --> Based on scree plot three dimensions will be sufficient for NMDS #

# MDS and plot stress using a Shepherd Plot #
MDS = metaMDS(GRIN, distance = "bray", k=3)
MDS$stress
stressplot(MDS) 
goodness(MDS)
# --> Shepherd plots showcase a not perfect, but acceptable R^2 value #

# Turn MDS points into a dataframe with treatment data for use in ggplot #
NMDS = data.frame(MDS = MDS$points, Treatment = Treat$Treatment, 
                  Plot = Treat$ID)

# NMDS Graphs
NMDS_plot = 
  ggplot() +
  geom_point(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, color = Treatment)) +
  #geom_text(data = species.scores, aes(x = NMDS1, y = NMDS2, label = species)) +
  stat_ellipse(data = NMDS, aes(x = MDS.MDS1, y = MDS.MDS2, color = Treatment), 
               linetype = "dashed", show.legend = T) +
  scale_color_manual(values=c("#FF3399", "#117733", "#3366FF")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family="Roboto Mono")) +
  labs(x = "MDS1", y = "MDS2", title = "NMDS Ordiantion - Seedbank Data")
NMDS_plot

ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/NMDS_Seedbank.png", width = 10, height = 7)

# Perform adonis to test the significance of treatments#

adon.results <- adonis2(GRIN ~ NMDS$Treat, method="bray",perm=999)
print(adon.results)

