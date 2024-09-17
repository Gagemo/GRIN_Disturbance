################################################################################
################################################################################
#########################      GRIN - Fire        ##############################
#########################    HEAT - Pres/Abs      ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
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
GRIN = read.csv("Data/GRIN - 2020-2023.csv")

GRIN$Coverage = as.numeric(GRIN$Coverage)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
GRIN <- mutate(GRIN, Coverage = case_when(
  grepl(10, Coverage) ~ 97.5,
  grepl(0, Coverage) ~ 0,
  grepl(1, Coverage) ~ 0.1,
  grepl(2, Coverage) ~ 0.5,
  grepl(3, Coverage) ~ 1.5,
  grepl(4, Coverage) ~ 3.5,
  grepl(5, Coverage) ~ 7.5,
  grepl(6, Coverage) ~ 17.5,
  grepl(7, Coverage) ~ 37.5,
  grepl(8, Coverage) ~ 62.5,
  grepl(9, Coverage) ~ 85
))

str(GRIN)
summary(GRIN)

# Selects just Seeding Treatment # 
GRIN = filter(GRIN, Treatment == 'S')

GRIN$YID <- paste(GRIN$Year,GRIN$ID)

GRIN$ID_ <- paste(GRIN$Fire, GRIN$ID)

#Renames values in fire treatments for heat map later #
GRIN$Fire <- recode(GRIN$Fire, 
                    C ="No Burn", W = "Winter", Sp = "Late-Spring")

# Separate Pre & Post-Treatment Data #
GRIN_22 = filter(GRIN, Year == "2")
GRIN_23 = filter(GRIN, Year == "3")

# Orders years and treatments so that they display in same sequence in graphs #
GRIN$Year = factor(GRIN$Year, levels=c('2','3'))
GRIN$Fire = factor(GRIN$Fire, levels=c('C','Sp','W'))

# Create Species Pivot Table with revisions to analyses #

# Select 2022 Data matrix #
Spp_22 <- dplyr::select(GRIN_22, ID_, Species, Coverage) %>% 
  matrify()

# Select 2022 Treatment Data matrix #
Spp_23 <- dplyr::select(GRIN_23, ID_, Species, Coverage) %>% 
  matrify()

# Creates treatment data frame for heat map #
Treat_22 = ungroup(GRIN_22) %>% 
  dplyr::select(ID_, Fire) %>%
  group_by(Fire, ID_) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID_')

Treat_23 = ungroup(GRIN_23) %>% 
  dplyr::select(ID_, Fire) %>%
  group_by(Fire, ID_) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID_')

ann_colors = list(
  Fire = c("No Burn" = "#333333", "Late-Spring" = "#FF9900", 
           "Winter" = "#3366FF"))


speciesHEAT_22 = pheatmap(Spp_22, show_rownames=F, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat_22, 
                       annotation_colors = ann_colors, fontsize = 20,
                       border_color = "black", display_numbers = FALSE,
                       cellheight=10, cellwidth = 20,
                       color=colorRampPalette(c("white", "orange", "red"))(50))
speciesHEAT_22

png(file = "Figures/Chapter 2 - Fire/PresAbs_22.png", 
    units="cm", width=30, height=30, res=100)
speciesHEAT_22
dev.off()

speciesHEAT_23 = pheatmap(Spp_23, show_rownames=F, cluster_cols=F, 
                          cluster_rows=F, annotation_row = Treat_23, 
                          annotation_colors = ann_colors, fontsize = 20,
                          border_color = "black", display_numbers = FALSE,
                          cellheight=10, cellwidth = 20,
                          color=colorRampPalette(c("white", "orange", "red"))(50))
speciesHEAT_23
png(file = "Figures/Chapter 2 - Fire/PresAbs_23.png", 
    units="cm", width=30, height=30, res=100)
speciesHEAT_23
dev.off()
