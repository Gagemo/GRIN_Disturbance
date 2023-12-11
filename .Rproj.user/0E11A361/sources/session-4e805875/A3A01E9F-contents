################################################################################
################################################################################
#########################   GRIN - Disturbance    ##############################
#########################    HEAT - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2021 - 2023         ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "labdsv", "pheatmap")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(pheatmap)

##########################     Read in  Data       #########################
GRIN = read.csv("Data/GRIN - 2021-2023.csv")

GRIN$Coverage = as.numeric(GRIN$Coverage)

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

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

GRIN = filter(GRIN, Species != "UK - White Forb")
GRIN = filter(GRIN, Species != "UK - Whirledleaf")
GRIN = filter(GRIN, Species != "UK - Mimosa")
GRIN = filter(GRIN, Species != "UK - Euphorb")

GRIN$YID <- paste(GRIN$Year,GRIN$ID)
GRIN$ID_ <- paste(GRIN$Treatment, GRIN$ID)

#Renames values in fire treatments for heat map later #
GRIN$Treatment <- recode(GRIN$Treatment, 
                         C ="No-Till", Tw = "Winter", Tsp = "Late-Spring")

# Seperate Pre & Post-Treatment Data #
GRIN_21 = filter(GRIN, Year == "1")
GRIN_22 = filter(GRIN, Year == "2")

# Orders years and treatments so that they display in same sequence in graphs #
GRIN$Year = factor(GRIN$Year, levels=c('1','2'))
GRIN$Treatment = factor(GRIN$Treatment, levels=c('C','Tsp','Tw'))

# Create Species Pivot Table with revisions to analyses #
# All treatments #
Spp <- dplyr::select(GRIN, ID_, Species, Coverage) %>% 
  matrify()

# Select 2021 Data matrix #
Spp_21 <- dplyr::select(GRIN_21, ID_, Species, Coverage) %>% 
  matrify()

# Select 2022 Treatment Data matrix #
Spp_22 <- dplyr::select(GRIN_22, ID_, Species, Coverage) %>% 
  matrify()

Treat = ungroup(GRIN_21) %>% 
  dplyr::select(ID_, Treatment) %>%
  group_by(Treatment, ID_) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID_')

ann_colors = list(
  Treatment = c("No-Till" = "#FF99CC", "Late-Spring" = "#FFFF99", 
                "Winter" = "#99CCFF"))

speciesHEAT = pheatmap(Spp_21, show_rownames=F, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat, 
                       annotation_colors = ann_colors, fontsize = 12,
                       border_color = "black", display_numbers = FALSE,
                       cellheight=10, cellwidth = 10,
                       color=colorRampPalette(c("white", "orange", "red"))(50))
speciesHEAT

png(file = "Figures/Chapter 1 - Soil Disturbance Seasonality/21_presAb.png",
    height = 2000, width = 1500, units="px", res = 150)
speciesHEAT
dev.off()

speciesHEAT = pheatmap(Spp_22, show_rownames=F, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat, 
                       annotation_colors = ann_colors, fontsize = 12,
                       border_color = "black", display_numbers = FALSE,
                       cellheight=10, cellwidth = 10,
                       color=colorRampPalette(c("white", "orange", "red"))(50))
speciesHEAT

png(file = "Figures/Chapter 1 - Soil Disturbance Seasonality/22_presAb.png",
    height = 2000, width = 1500, units="px", res = 150)
speciesHEAT
dev.off()

