######################  GRIN - Soil Disturbance Seasonality ####################
######################        Species Richness              ####################
######################      University of Florida           ####################
######################          Gage LaPierre               ####################
######################             2022                     ####################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "labdsv", "reshape")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################

library(tidyverse)
library(vegan)
library(labdsv)
library(reshape)

##########################     Read in 2022 Data       ##############################
data = read.csv("Data/GRIN - 2022.csv")
data$Coverage = as.numeric(data$Coverage)

str(data)
summary(data)

### Remove NA and empty Values ###
data = filter(data, Coverage != "NA") %>%
  filter(Coverage != "")

# Remove Seeding Treatment # 
data = filter(data, Treatment != 'S')

# Create Species Pivot Table #
Spp = dplyr::select(data, ID, Species, Coverage) %>% matrify() 
Spp[] <- lapply(Spp, as.numeric)

# Create Grouped Treatment/ Environment Table and Summaries to fit Species Table #
Treat <- select(data, Treatment, ID)%>% group_by(ID, Treatment) %>% summarise()

SR = specnumber(Spp)
SR = as.data.frame(SR)

## Merge species richness with habitat/plot data for ggplot ##
SR_treat = cbind(Treat, SR) 

## Species Richness Boxplot ##
SR_Box = 
  ggplot(SR_treat, aes(x = Treatment, y = SR, fill  = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", alpha=0.7, width = 0.25) +
  scale_fill_manual(values=c("#FF3399", "#FFFF33", "#3366FF"))+
  labs(x="", y = "Species Richness") +
  theme_classic() +
  theme(legend.position = "none")
SR_Box
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/SR_Box_2022.png")


SR_anova = aov(SR ~ Treatment, data = SR_treat)
summary(SR_anova)
tukey.one.way<-TukeyHSD(SR_anova)
tukey.one.way

##########################     Read in 2021 Data       ##############################
data = read.csv("Data/GRIN - 2021.csv")
data$Coverage = as.numeric(data$Coverage)

str(data)
summary(data)

### Remove NA and empty Values ###
data = filter(data, Coverage != "NA") %>%
  filter(Coverage != "")

# Remove Seeding Treatment # 
data = filter(data, Treatment != 'S')

# Create Species Pivot Table #
Spp = dplyr::select(data, ID, Species, Coverage) %>% matrify() 
Spp[] <- lapply(Spp, as.numeric)

# Create Grouped Treatment/ Environment Table and Summaries to fit Species Table #
Treat <- select(data, Treatment, ID)%>% group_by(ID, Treatment) %>% summarise()

SR = specnumber(Spp)
SR = as.data.frame(SR)

## Merge species richness with habitat/plot data for ggplot ##
SR_treat = cbind(Treat, SR) 

## Species Richness Boxplot ##
SR_Box = 
  ggplot(SR_treat, aes(x = Treatment, y = SR, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", alpha=0.7, width = 0.25)  +
  scale_fill_manual(values=c("#FF3399", "#FFFF33", "#3366FF")) +
  labs(x="", y = "Species Richness") +
  theme_classic() +
  theme(legend.position = "none")
SR_Box
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/SR_Box_2021.png")


SR_anova = aov(SR ~ Treatment, data = SR_treat)
summary(SR_anova)
tukey.one.way<-TukeyHSD(SR_anova)
tukey.one.way
