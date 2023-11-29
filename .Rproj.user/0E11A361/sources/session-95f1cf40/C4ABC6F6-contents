#########################   GRIN - Disturbance    ##############################
#########################       Bahaiagrass       ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
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

##########################Read in 2021 - 2023 Data  ############################

GRIN = read.csv("Data/GRIN - 2021-2023.csv")
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

# Creates data sets by year #
PN_21 = filter(PN, Year == 1)
PN_22 = filter(PN, Year == 2)
PN_23 = filter(PN, Year == 3)

## Bareground Coverage ##
box = 
  ggplot(PN, aes(x = Treatment, y = Coverage, fill = Treatment)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(vars(Year)) +
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

ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2021-2023_BahaiaGrass.png", 
       width = 10, height = 7)

# Test for Significance across years #
anova = aov(Coverage ~ Treatment, data = PN_21)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

anova = aov(Coverage ~ Treatment, data = PN_22)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

anova = aov(Coverage ~ Treatment, data = PN_23)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

## Bareground Coverage ##
box = 
  ggplot(PN, aes(x = Treatment, y = Coverage, fill = Treatment)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(size=3, alpha = 0.5, color="black", width = 0.25) +
  geom_signif(y_position = c(100, 100), xmin = c(0.7,2), xmax = c(1.9,3),
              annotation=c("***", "***"), size = 0.8,
              textsize = 5, tip_length = 0.01) +
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