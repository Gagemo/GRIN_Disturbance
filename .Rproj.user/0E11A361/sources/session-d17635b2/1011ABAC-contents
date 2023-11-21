#########################   GRIN - Disturbance    ##############################
#########################    HEAT - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2021 - 2023         ##############################
################################################################################
################################################################################
################################################################################
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
GRIN = filter(GRIN, Treatment != 'Tw')
GRIN = filter(GRIN, Treatment != 'Tsp')

GRIN$YID <- paste(GRIN$Year,GRIN$ID)

# Seperate Pre & Post-Treatment Data #
GRIN_22 = filter(GRIN, Year == "2")
GRIN_23 = filter(GRIN, Year == "3")

# Orders years and treatments so that they display in same sequence in graphs #
GRIN$Year = factor(GRIN$Year, levels=c('2','3'))
GRIN$Fire = factor(GRIN$Fire, levels=c('C','Sp','W'))
GRIN$Treatment = factor(GRIN$Treatment, levels=c('C', 'S'))

# Create Species Pivot Table with Rae's revisions to analyses #
# All treatments #
Spp <- dplyr::select(GRIN, YID, Species, Coverage) %>% 
  matrify()

# Select 2021 Data matrix #
Spp_22 <- dplyr::select(GRIN_22, YID, Species, Coverage) %>% 
  matrify()

# Select 2022 Treatment Data matrix #
Spp_23 <- dplyr::select(GRIN_23, YID, Species, Coverage) %>% 
  matrify()

#################### Species abundances ########################################
# Creates and joins  data year 22 & 23 to make long data format #
Two_Abundance <- GRIN[which(GRIN$Year == "2"),]
Three_Abundance <- GRIN[which(GRIN$Year == "3"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance, 
                         by = c('ID', 'Treatment', "Fire", 'Species'))

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundace to reflect percentage change from pretreatment (Year 2) to Post Treatment (Year 3)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(ID, Treatment, Fire, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(ID, Treatment, Fire, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x)

Treat = ungroup(Change_Abundance) %>% 
  dplyr::select(ID, Treatment, Fire) %>%
  group_by(ID, Treatment, Fire) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID')

GRIN_change <- ungroup(Change_Abundance) %>%
  dplyr::select(ID, Species, Change_abundance) %>%
  as.data.frame() %>%
  matrify()

GRIN_change = as.matrix(GRIN_change[, -1])

ann_colors = list(
  Fire = c("C" = "#FF3399", "Sp" = "#117733", "W" = "#3366FF"),
  Treatment = c("C"  = "grey", "S" = "yellow"))

speciesHEAT = pheatmap(GRIN_change, show_rownames=F, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat, cex = 1, 
                       annotation_colors = ann_colors, fontsize = 12,
)
speciesHEAT
png(file = "Figures/Chapter 2 - Fire/Heat.png", 
    units="cm", width=40, height=22, res=999)
speciesHEAT
dev.off()