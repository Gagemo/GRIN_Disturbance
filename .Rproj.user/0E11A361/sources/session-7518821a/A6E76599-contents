################################################################################
################################################################################
#########################       GRIN - Fire       ##############################
#########################      Cover Change       ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2020 - 2023        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggsignif", "multcompView", "ggpubr", "rstatix",
                      "vegan", "labdsv")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(extrafont)
#font_import()
loadfonts(device = "win")
library(tidyverse)
library(vegan)
library(agricolae)
library(ggsignif)
library(multcompView)
library(ggpubr)
library(rstatix)
library(vegan)
library(labdsv)

##########################     Read in 2020 - 2023 Data       ##################
GRIN = read.csv("Data/GRIN_FUN-2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Select just Seeding Treatment # 
GRIN = filter(GRIN, Treatment == 'S')

# Selects just seeded species #
GRIN = filter(GRIN, Species == "Aristida stricta" |
                Species == "Andropogon ternarius" |
                Species == "Eragrostis spectabilis" |
                Species == "Liatris gracilis" |
                Species == "Pityopsis graminifolia" |
                Species == "Sorgastrum secundum")

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

GRIN$Coverage = as.numeric(GRIN$Coverage)

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

# Creates data sets by year #
GRIN_22 = filter(GRIN, Year == 2)
GRIN_23 = filter(GRIN, Year == 3)

################################################################################
########################## Seeded Species Analysis #############################
################################################################################

# Creates and joins  data year 22 & 23 to make long data format #
Two_Abundance <- GRIN[which(GRIN$Year == "2"),]
Three_Abundance <- GRIN[which(GRIN$Year == "3"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance, 
                         by = c('ID', "Fire", 'Species'))
Abundance_w = arrange(Abundance_w, Fire)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(ID, Fire, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(ID, Fire, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x) %>%
  filter(Change_abundance != 0)

change = 
  ggplot(Change_Abundance, aes(x = Fire, y = Change_abundance), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  facet_wrap(vars(Species)) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Total (%) Coverage", title = "")
change
ggsave("Figures/Chapter 2 - Fire/change.png", 
       width = 10, height = 7)

seed = 
  ggplot(GRIN, aes(x = Fire, y = Coverage)) +
  geom_boxplot(aes(fill=Species), alpha = 0.5, outlier.shape = NA) +
  facet_wrap(vars(Year)) +
  geom_point(aes(fill=Species), position = position_jitterdodge(), alpha = 0.5) +
  #scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
  #                  values=c("#333333", "#FF9900", "#3366FF")) +
  #scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "bottom") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Total (%) Coverage", title = "")
seed
ggsave("Figures/Chapter 2 - Fire/seed.png", 
       width = 10, height = 7)
