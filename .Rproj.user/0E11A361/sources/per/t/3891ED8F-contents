################################################################################
################################################################################
#########################      GRIN - Fire        ##############################
#########################       Lovegrass         ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
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
GRIN = filter(GRIN, Treatment != 'Tsp')
GRIN = filter(GRIN, Treatment != 'Tw')
GRIN = filter(GRIN, Treatment != 'C')

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

GRIN = filter(GRIN, Species == "Eragrostis spectabilis")
summary(GRIN)

#Renames values in fire treatments for heat map later #
GRIN$Fire <- recode(GRIN$Fire, Sp="Late-Spring", W = "Winter", C = "No Burn")

# Creates data sets by year #
GRIN_21 = filter(GRIN, Year == 1)
GRIN_22 = filter(GRIN, Year == 2)
GRIN_23 = filter(GRIN, Year == 3)

# Label Years for ggplot #
year_names <- as_labeller(
  c(`1` = "2021", `2` = "2022", `3` = "2023"))

## Bareground Coverage ##
box = 
  ggplot(GRIN, aes(x = Fire, y = Coverage, fill = Fire)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(.~Year, labeller = year_names) +
  geom_jitter(size=3, alpha = 0.5, color="black", width = 0.25) +
  scale_fill_manual(values=c("#333333", "#FF9900", "#3366FF")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family = "Roboto Mono"))+
  theme_classic() 
box

ggsave("Figures/Chapter 2 - Fire/2021-2023_Lovegrass.png", 
       width = 10, height = 7)

# Test for Significance across years #
anova = aov(Coverage ~ Fire, data = GRIN_21)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

anova = aov(Coverage ~ Fire, data = GRIN_22)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way

anova = aov(Coverage ~ Fire, data = GRIN_23)
summary(anova)
tukey.one.way<-TukeyHSD(anova)
tukey.one.way
