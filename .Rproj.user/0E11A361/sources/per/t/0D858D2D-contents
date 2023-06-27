#########################   GRIN - Disturbance    ##############################
#########################       Bareground        ##############################
#########################  University of Florida  ##############################
#########################       Gage LaPierre     ##############################
#########################          2021 - 2022    ##############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "agricolae")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################

library(tidyverse)
library(vegan)
library(agricolae)

##########################     Read in 2022 Data       ##############################

data = read.csv("Data/GRIN - 2022.csv")
data$Coverage = as.numeric(data$Coverage)
data$Plot = as.character(data$Plot)

str(data)
summary(data)

# Remove Seeding Treatment # 
data = filter(data, Treatment != 'S')

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
data <- mutate(data, Coverage = case_when(
  grepl(1, Coverage) ~ 0.1,
  grepl(2, Coverage) ~ 0.5,
  grepl(3, Coverage) ~ 1.5,
  grepl(4, Coverage) ~ 3.5,
  grepl(5, Coverage) ~ 7.5,
  grepl(6, Coverage) ~ 17.5,
  grepl(7, Coverage) ~ 37.5,
  grepl(8, Coverage) ~ 62.5,
  grepl(9, Coverage) ~ 85,
  grepl(10, Coverage) ~ 97.5
))

bare = filter(data, Group == "Bare")

## Bareground Coverage ##
box = 
  ggplot(bare, aes(x = Treatment, y = Coverage, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", alpha=0.7, width = 0.25) +
  scale_fill_manual(values=c("#FF3399", "#FFFF33", "#3366FF"))+
  theme_classic() 
box

ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2022_Bareground.png", width = 10, height = 7)

# Test for Significance #
anova = aov(Coverage ~ Treatment, data = data)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

##########################     Read in 2021 Data       #########################
data = read.csv("Data/GRIN - 2021.csv")
data$Coverage = as.numeric(data$Coverage)
data$Plot = as.character(data$Plot)

str(data)
summary(data)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
data <- mutate(data, Coverage = case_when(
  grepl(1, Coverage) ~ 0.1,
  grepl(2, Coverage) ~ 0.5,
  grepl(3, Coverage) ~ 1.5,
  grepl(4, Coverage) ~ 3.5,
  grepl(5, Coverage) ~ 7.5,
  grepl(6, Coverage) ~ 17.5,
  grepl(7, Coverage) ~ 37.5,
  grepl(8, Coverage) ~ 62.5,
  grepl(9, Coverage) ~ 85,
  grepl(10, Coverage) ~ 97.5
))

bare = filter(GRIN, Group == "Bare")

## Bareground Coverage ##

box = 
  ggplot(bare, aes(x = Treatment, y = Coverage, fill = Treatment)) +
  geom_boxplot() +
  geom_jitter(color="black", alpha=0.7, width = 0.25) +
  scale_fill_manual(values=c("#FF3399", "#66FF33", "#FFFF33", "#3366FF"))+
  theme_classic() 
box

ggsave("Figures/2021_Bareground.png", width = 10, height = 7)

# Test for Significance #
anova = aov(Coverage ~ Treatment, data = GRIN)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way
