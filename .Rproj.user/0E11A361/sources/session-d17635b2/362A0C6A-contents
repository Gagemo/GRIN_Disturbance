#########################   GRIN - Disturbance    ##############################
#########################    Functional Groups    ##############################
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

#########################     Installs Packages   ##############################

list.of.packages <- c("tidyverse", "vegan", "agricolae")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################

library(tidyverse)
library(vegan)
library(agricolae)

##########################     Read in 2021-2023 Data       ####################
GRIN = read.csv("Data/GRIN - 2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

GRIN <- GRIN %>% mutate(Life = case_when(
    grepl("Andropogon ternarius", Species) ~ "Perennial",
    grepl("Andropogon virginicus", Species) ~ "Perennial",
    grepl("Aristida stricta", Species) ~ "Perennial",
    grepl("Baccharis halimifolia", Species) ~ "Perennial",
    grepl("Bulbostylis spp.", Species) ~ "Annual",
    grepl("Cirsium spp.", Species) ~ "Perennial",
    grepl("Crotalaria rotundifolia", Species) ~ "Perennial",
    grepl("Crotalaria spectabilis", Species) ~ "Perennial",
    grepl("Cynodon dactylon", Species) ~ "Perennial",
    grepl("Cyperus esculentus", Species) ~ "Perennial",
    grepl("Cyperus spp.", Species) ~ "Perennial",
    grepl("Dichanthelium spp.", Species) ~ "Perennial",
    grepl("Dichondra spp.", Species) ~ "Perennial",
    grepl("Digitaria spp.", Species) ~ "Annual",
    grepl("Eragrostis spectabilis", Species) ~ "Perennial",
    grepl("Eremochloa ophiuroides", Species) ~ "Perennial",
    grepl("Erigeron canadensis", Species) ~ "Annual",
    grepl("Eupatorium capilifolium", Species) ~ "Perennial",
    grepl("Galactia spp.", Species) ~ "Perennial",
    grepl("Gamochaeta spp.", Species) ~ "Annual",
    grepl("Indigofera hirsuta", Species) ~ "Perennial",
    grepl("Lactuca canadensis", Species) ~ "Perennial",
    grepl("Liatris gracilis", Species) ~ "Perennial",
    grepl("Liquidambar styraciflua", Species) ~ "Perennial",
    grepl("Oenothera laciniata", Species) ~ "Perennial",
    grepl("Oxalis corniculata", Species) ~ "Annual",
    grepl("Passiflora spp.", Species) ~ "Perennial",
    grepl("Paspalum notatum", Species) ~ "Perennial",
    grepl("Paspalum setaceum", Species) ~ "Perennial",
    grepl("Phyllanthus urinaria", Species) ~ "Annual",
    grepl("Phytolacca americana", Species) ~ "Perennial",
    grepl("Pityopsis graminifolia", Species) ~ "Perennial",
    grepl("Richardia spp.", Species) ~ "Annual",
    grepl("Rubus spp.", Species) ~ "Perennial",
    grepl("Rumex hastatulus", Species) ~ "Annual",
    grepl("Solanum viarum", Species) ~ "Perennial",
    grepl("Sorghastrum secundum", Species) ~ "Perennial",
    grepl("Polypremum procumbens", Species) ~ "Perennial",
    grepl("UK - Blue Vervain", Species) ~ "Perennial",
    grepl("Dalea pinnata", Species) ~ "Perennial",
    grepl("Raphanus raphanistrum", Species) ~ "Annual",
    grepl("Urtica dioica", Species) ~ "Perennial",
    grepl("Verbena spp.", Species) ~ "Annual",
    TRUE ~ ""
  ))

# Plant Family #

GRIN <- GRIN %>% 
  mutate(Family = case_when(
    grepl("Andropogon ternarius", Species) ~ "Poaceae",
    grepl("Andropogon virginicus", Species) ~ "Poaceae",
    grepl("Aristida stricta", Species) ~ "Poaceae",
    grepl("Baccharis halimifolia", Species) ~ "Asteraceae",
    grepl("Bulbostylis spp.", Species) ~ "Cyperaceae",
    grepl("Cirsium spp.", Species) ~ "Asteraceae",
    grepl("Crotalaria rotundifolia", Species) ~ "Fabaceae",
    grepl("Crotalaria spectabilis", Species) ~ "Fabaceae",
    grepl("Cynodon dactylon", Species) ~ "Poaceae",
    grepl("Cyperus esculentus", Species) ~ "Cyperaceae",
    grepl("Cyperus spp.", Species) ~ "Cyperaceae",
    grepl("Dichanthelium spp.", Species) ~ "Poaceae",
    grepl("Dichondra spp.", Species) ~ "Convolvulaceae",
    grepl("Digitaria spp.", Species) ~ "Poaceae",
    grepl("Eragrostis spectabilis", Species) ~ "Poaceae",
    grepl("Eremochloa ophiuroides", Species) ~ "Poaceae",
    grepl("Erigeron canadensis", Species) ~ "Asteraceae",
    grepl("Eupatorium capilifolium", Species) ~ "Asteraceae",
    grepl("Galactia spp.", Species) ~ "Fabaceae",
    grepl("Gamochaeta spp.", Species) ~ "Asteraceae",
    grepl("Indigofera hirsuta", Species) ~ "Fabaceae",
    grepl("Lactuca canadensis", Species) ~ "Asteraceae",
    grepl("Liatris gracilis", Species) ~ "Asteraceae",
    grepl("Liquidambar styraciflua", Species) ~ "Altingiaceae",
    grepl("Oenothera laciniata", Species) ~ "Onagraceae",
    grepl("Oxalis corniculata", Species) ~ "Oxalidaceae",
    grepl("Passiflora spp.", Species) ~ "Passifloraceae",
    grepl("Paspalum notatum", Species) ~ "Poaceae",
    grepl("Paspalum setaceum", Species) ~ "Poaceae",
    grepl("Phyllanthus urinaria", Species) ~ "Fabaceae",
    grepl("Phytolacca americana", Species) ~ "Phytolaccaceae",
    grepl("Pityopsis graminifolia", Species) ~ "Asteraceae",
    grepl("Richardia spp.", Species) ~ "Rubiaceae",
    grepl("Rubus spp.", Species) ~ "Asteraceae",
    grepl("Rumex hastatulus", Species) ~ "Polygonaceae",
    grepl("Solanum viarum", Species) ~ "Solanaceae",
    grepl("Sorghastrum secundum", Species) ~ "Poaceae",
    grepl("Polypremum procumbens", Species) ~ "Tetrachondraceae",
    grepl("UK - Blue Vervain", Species) ~ "Asteraceae",
    grepl("Dalea pinnata", Species) ~ "Asteraceae",
    grepl("Raphanus raphanistrum", Species) ~ "Brassicaceae",
    grepl("Urtica dioica", Species) ~ "Urticaceae",
    grepl("Verbena spp.", Species) ~ "	Verbanaceae",
    TRUE ~ ""
  ))

# Raunkier #

GRIN <- GRIN %>%  
  mutate(Raunkier = case_when(
    grepl("Andropogon ternarius", Species) ~ "Hemicryptophyte",
    grepl("Andropogon virginicus", Species) ~ "Hemicryptophyte",
    grepl("Aristida stricta", Species) ~ "Hemicryptophyte",
    grepl("Baccharis halimifolia", Species) ~ "Phanerophyte",
    grepl("Bulbostylis spp.", Species) ~ "Therophyte",
    grepl("Cirsium spp.", Species) ~ "Hemicryptophyte",
    grepl("Crotalaria rotundifolia", Species) ~ "Hemicryptophyte",
    grepl("Crotalaria spectabilis", Species) ~ "Hemicryptophyte",
    grepl("Cynodon dactylon", Species) ~ "Hemicryptophyte",
    grepl("Cyperus esculentus", Species) ~ "Geophyte",
    grepl("Cyperus spp.", Species) ~ "Hemicryptophyte",
    grepl("Dichanthelium spp.", Species) ~ "Hemicryptophyte",
    grepl("Dichondra spp.", Species) ~ "Hemicryptophyte",
    grepl("Digitaria spp.", Species) ~ "Therophyte",
    grepl("Eragrostis spectabilis", Species) ~ "Hemicryptophyte",
    grepl("Eremochloa ophiuroides", Species) ~ "Hemicryptophyte",
    grepl("Erigeron canadensis", Species) ~ "Therophyte",
    grepl("Eupatorium capilifolium", Species) ~ "Hemicryptophyte",
    grepl("Galactia spp.", Species) ~ "Hemicryptophyte",
    grepl("Gamochaeta spp.", Species) ~ "Therophyte",
    grepl("Indigofera hirsuta", Species) ~ "Hemicryptophyte",
    grepl("Lactuca canadensis", Species) ~ "Therophyte",
    grepl("Liatris gracilis", Species) ~ "Geophyte",
    grepl("Liquidambar styraciflua", Species) ~ "Phanerophyte",
    grepl("Oenothera laciniata", Species) ~ "Hemicryptophyte",
    grepl("Oxalis corniculata", Species) ~ "Therophyte",
    grepl("Passiflora spp.", Species) ~ "Geophyte",
    grepl("Paspalum notatum", Species) ~ "Hemicryptophyte",
    grepl("Paspalum setaceum", Species) ~ "Hemicryptophyte",
    grepl("Phyllanthus urinaria", Species) ~ "Therophyte",
    grepl("Phytolacca americana", Species) ~ "Hemicryptophyte",
    grepl("Pityopsis graminifolia", Species) ~ "Hemicryptophyte",
    grepl("Richardia spp.", Species) ~ "Therophyte",
    grepl("Rubus spp.", Species) ~ "Chamaephyte",
    grepl("Rumex hastatulus", Species) ~ "Therophyte",
    grepl("Solanum viarum", Species) ~ "Hemicryptophyte",
    grepl("Sorghastrum secundum", Species) ~ "Hemicryptophyte",
    grepl("Polypremum procumbens", Species) ~ "Hemicryptophyte",
    grepl("UK - Blue Vervain", Species) ~ "Hemicryptophyte",
    grepl("Dalea pinnata", Species) ~ "Hemicryptophyte",
    grepl("Raphanus raphanistrum", Species) ~ "Therophyte",
    grepl("Urtica dioica", Species) ~ "Hemicryptophyte",
    grepl("Verbena spp.", Species) ~ "Therophyte",
    TRUE ~ ""
  ))

write.csv(GRIN, "Data/GRIN_FUN-2021-2023.csv", row.names=FALSE)
