#########################   GRIN - Disturbance    ##############################
#########################    Species Richness     ##############################
#########################  University of Florida  ##############################
#########################       Gage LaPierre     ##############################
#########################          2022           ##############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################

library(tidyverse)
library(vegan)

##########################     Read in Data       ##############################
GRIN = read.csv("Data/GRIN.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)

str(GRIN)
summary(GRIN)

# Create Species Pivot Table #
Spp = dplyr::select(GRIN, ID, Species, Coverage) %>% matrify() 
Spp[] <- lapply(Spp, as.numeric)

# Create grouped table & summaries to fit species table #
Treat = group_by(GRIN, ID, Treatment) %>% summarise()

# Species Richness #
table_SR <- table(GRIN$Species, GRIN$ID)
table_SR 

SR = specnumber(table_SR , MARGIN=2)
SR = as.data.frame(SR)

## Merge species richness with habitat/plot data for ggplot ##
SR_treat = cbind(Treat, SR) 

## Species Richness Boxplot ##
SR_Box = 
  ggplot(SR_treat, aes(x = Treatment, y = SR, color  = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", alpha=0.7, width = 0.25) +
  labs(x="", y = "Species Richness") +
  theme_classic() +
  theme(legend.position = "none")
SR_Box
ggsave("Figures/SR_Box.png")


SR_anova = aov(SR ~ Treatment, data = SR_treat)
summary(SR_anova)
tukey.one.way<-TukeyHSD(SR_anova)
tukey.one.way
