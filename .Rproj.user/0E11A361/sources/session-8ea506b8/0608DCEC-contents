#########################   GRIN - Disturbance    ##############################
#########################    NMDS - Community     ##############################
######################### University of Florida   ##############################
#########################    Gage LaPierre        ##############################
#########################     2021 - 2022         ##############################
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
GRIN = read.csv("Data/GRIN - 21-22.csv")

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

GRIN$YID <- paste(GRIN$Year,GRIN$ID)

# Seperate Pre & Post-Treatment Data #
GRIN_21 = filter(GRIN, Year == "2021")
GRIN_22 = filter(GRIN, Year == "2022")

# Orders years and treatments so that they display in same sequence in graphs #
GRIN$Year = factor(GRIN$Year, levels=c('2021','2022'))
GRIN$Fire = factor(GRIN$Treatment, levels=c('C','Tsp','Tw'))

# Create Species Pivot Table with Rae's revisions to analyses #
# All treatments #
Spp <- dplyr::select(GRIN, YID, Species, Coverage) %>% 
  matrify()

# Select 2021 Data matrix #
Spp_21 <- dplyr::select(GRIN_21, YID, Species, Coverage) %>% 
  matrify()

# Select 2022 Treatment Data matrix #
Spp_22 <- dplyr::select(GRIN_22, YID, Species, Coverage) %>% 
  matrify()

#################### Species abundances ########################################
# Creates and joins  data year 21 & 22 to make long data format #
One_Abundance <- GRIN[which(GRIN$Year == "2021"),]
Three_Abundance <- GRIN[which(GRIN$Year == "2022"),]

Abundance_w <- full_join(One_Abundance, Three_Abundance, 
                         by = c('ID', 'Treatment', 'Species'))

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                    Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                    Abundance_w$Coverage.y)

# Change abundace to reflect percentage change from pretreatment (Year 1) to Post Treatment (Year 3)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(ID, Treatment, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(ID, Treatment, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x)

Treat = ungroup(Change_Abundance) %>% 
  dplyr::select(ID, Treatment) %>%
  group_by(ID, Treatment) %>%
  summarise() %>%
  remove_rownames() %>%
  column_to_rownames(var = 'ID')

GRIN_change <- ungroup(Change_Abundance) %>%
  dplyr::select(ID, Species, Change_abundance) %>%
  as.data.frame() %>%
  matrify()

GRIN_change = as.matrix(GRIN_change[, -1])

ann_colors = list(
  Treatment = c("C" = "#FF3399", "Tsp" = "#117733", "Tw" = "#3366FF"))

speciesHEAT = pheatmap(GRIN_change, show_rownames=F, cluster_cols=F, 
                       cluster_rows=F, annotation_row = Treat, cex = 1, 
                       annotation_colors = ann_colors, fontsize = 12,
                       )
speciesHEAT
png(file = "Figures/Chapter 1 - Soil Disturbance Seasonality/Heat.png", 
    units="cm", width=40, height=22, res=999)
speciesHEAT
dev.off()