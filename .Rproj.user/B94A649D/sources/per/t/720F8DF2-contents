#########################   GRIN - Disturbance    ##############################
#########################    Functional Groups    ##############################
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

list.of.packages <- c("tidyverse", "vegan", "agricolae")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################

library(tidyverse)
library(vegan)
library(agricolae)

##########################     Read in Data       ##############################
GRIN = read.csv("Data/GRIN.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
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
## Average Coverage by Functional Group ##

Fun = GRIN %>% group_by(Group, ID, Plot, Sub_Plot, Treatment) %>% 
  summarise(Coverage = sum(Coverage))

## Coverage by Functional Group per Plot Boxplot ##
Fun_Box = 
  ggplot(Fun, aes(x = Group, y = Coverage, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.25)+
  facet_wrap(vars(Treatment)) +
  theme_classic() 
Fun_Box

# Test for Significance #
Fun_anova = aov(Coverage ~ Treatment + Group, data = Fun)
summary(Fun_anova)
tukey.one.way<-TukeyHSD(Fun_anova)
tukey.one.way
HSD.stat = HSD.test(Fun_anova,trt = c("Treatment", "Group")) #HSD Tukey 
HSD.stat

