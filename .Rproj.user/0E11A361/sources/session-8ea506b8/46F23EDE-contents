#########################   GRIN - Disturbance    ##############################
#########################    Functional Groups    ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2022 - 2023        ##############################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggsignif")
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

##########################     Read in 2022 Data       #########################
GRIN = read.csv("Data/GRIN_FUN-2022.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

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

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

# Composition of Treatments using mean coverage per species. 
df = group_by(GRIN, Treatment, Group, Species) %>% 
  dplyr::summarize(Avg=mean(Coverage))
df = filter(df, Avg > 5) # Species <5% Coverage not included. 

Veg_Bar = 
  ggplot(df, aes(x = Treatment, y = Avg, fill = Group)) +
  geom_bar(position = "fill", stat = "identity", color = "black", alpha = 0.5) +
  geom_text(aes(label = Species), stat = "identity", 
            size = 5, position=position_fill(0.5), colour = "black",
            family = "Roboto Mono") +
  scale_fill_manual(breaks=c('Bare', 'Forb', 'Grass', 'Sedge'),
                    values = c("#993300", "#660066", "#339966", "#E7B800"), 
                    drop = FALSE, 
                    labels=c("Bare", "Forb", "Grass", "Sedge")) +
  geom_signif(y_position = c(1.01,1.01), xmin = c(0.7,2), xmax = c(1.9,3),
              annotation=c("***", "***"), tip_length = 0.0001) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family = "Roboto Mono"))+
  labs(x = "Treatment", y = "Coverage")
Veg_Bar
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2022_CompBar.png", 
       width = 10, height = 7)

################################################################################
################################################################################
################################ 2021 Data #####################################
################################################################################
################################################################################

##########################     Read in 2021 Data       #########################
GRIN = read.csv("Data/GRIN_FUN-2021.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

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

# Composition of Treatments using mean coverage per species. 
df = group_by(GRIN, Treatment, Group, Species) %>% 
  dplyr::summarize(Avg=mean(Coverage))
df = filter(df, Avg > 5) # Species <5% Coverage not included. 

Veg_Bar = 
  ggplot(df, aes(x = Treatment, y = Avg, fill = Group)) +
  geom_bar(position = "fill", stat = "identity", color = "black", alpha = 0.5) +
  geom_text(aes(label = Species), stat = "identity", 
            size = 4, position=position_fill(0.5), colour = "black",
            family = "Roboto Mono") +
  scale_fill_manual(values = c("#993300", "#660066", "#339966", "#E7B800")) +
  geom_signif(y_position = c(1.01,1.01), xmin = c(0.7,2), xmax = c(1.9,3),
              annotation=c("***", "***"), tip_length = 0.0001) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16,  family = "Roboto Mono")) +
  labs(x = "Treatment", y = "Coverage")
Veg_Bar
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2021_CompBar.png", 
       width = 10, height = 7)