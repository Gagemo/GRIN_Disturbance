#########################   GRIN - Disturbance    ##############################
#########################       Bareground        ##############################
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
#
#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "agricolae")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################

library(tidyverse)
library(vegan)
library(agricolae)

##########################     Read in 2022 Data  ##############################

GRIN = read.csv("Data/GRIN - 2022.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
  grepl(0, Coverage) ~ 0,
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

PN = filter(GRIN, Species == "Paspalum notatum")
summary(PN)
## Bareground Coverage ##
box = 
  ggplot(PN, aes(x = Treatment, y = Coverage, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size=3, alpha = 0.5, color="black", width = 0.25) +
  scale_fill_manual(values=c("#FF3399", "#117733", "#3366FF")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family = "Roboto Mono"))+
  theme_classic() 
box

ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2022_BahaiaGrass.png", 
       width = 10, height = 7)

# Test for Significance #
anova = aov(Coverage ~ Treatment, data = PN)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

##########################     Read in 2021 Data       #########################
GRIN = read.csv("Data/GRIN - 2021.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
  grepl(0, Coverage) ~ 0,
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

PN = filter(GRIN, Species == "Paspalum notatum")

## Paspalum notatum Coverage ##

box = 
  ggplot(PN, aes(x = Treatment, y = Coverage, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  geom_jitter(size=3, alpha = 0.5, color="black", width = 0.25) +
  scale_fill_manual(values=c("#FF3399", "#117733", "#3366FF")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family = "Roboto Mono"))+
  theme_classic() 
box
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2021_BahaiaGrass.png", 
       width = 10, height = 7)

# Test for Significance #
anova = aov(Coverage ~ Treatment, data = PN)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way
