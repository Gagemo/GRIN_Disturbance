################################################################################
################################################################################
#########################   GRIN - Disturbance    ##############################
#########################    Functional Groups    ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2022 - 2023        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggsignif", "multcompView")
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

##########################     Read in 2020 - 2023 Data       ##################
GRIN = read.csv("Data/GRIN_FUN-2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Year != 3)

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

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

# Sums CV data #
GRIN1 <- GRIN %>%
  select(ID, Year, Treatment, Group, Species, Coverage) %>% 
  group_by(ID, Year, Treatment, Group)  %>% 
  dplyr::summarize(TotalCV = sum(Coverage))

# Creates data sets by year #
GRIN_21 = filter(GRIN1, Year == 1)
GRIN_22 = filter(GRIN1, Year == 2)

################################################################################
##################### Functional Groups Analysis & Plot ########################
################################################################################
########################### 2021 Data Set ######################################
################################################################################

# Coverage by Functional Group per Plot Box Plot #
Group_Box = 
  ggplot(GRIN1, aes(x = Treatment, y = TotalCV, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA ) +
  facet_wrap(.~Year) +
  scale_fill_manual(breaks=c('Bare', 'Forb', 'Grass', 'Sedge', 'Woody'),
                    values = c("tan", "#660066", "#339966", "#E7B800", "brown"), 
                    labels=c("Bare", "Forb", "Grass", "Sedge", 'Woody')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="bold", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold")) +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "Treatment", y = "Total (%) Coverage")
Group_Box

###############################  Forbs 2021 ####################################
Forb = filter(GRIN_21, Group == "Forb")

# Test for Significance #
Fun_anova = aov(TotalCV ~ Treatment, data = Forb)
summary(Fun_anova)

tukey.one.way<-TukeyHSD(Fun_anova)
tukey.one.way

###############################  Forbs 2022 ####################################
Forb = filter(GRIN_22, Group == "Forb")

# Test for Significance #
Fun_anova = aov(TotalCV ~ Treatment, data = Forb)
summary(Fun_anova)

tukey.one.way<-TukeyHSD(Fun_anova)
tukey.one.way

################################################################################
############################## Annual vs Perennial #############################
################################################################################

LIFE <- GRIN %>% 
  filter(Life != "") %>%
  select(ID, Treatment, Life, Species, Coverage) %>% 
  group_by(ID, Treatment, Life)  %>% 
  summarize(TotalCV = sum(Coverage))

# Test for Significance #
Life_anova = aov(TotalCV ~ Treatment + Life, data = LIFE)
summary(Life_anova)
tukey.one.way<-TukeyHSD(Life_anova)
tukey.one.way
HSD.stat = HSD.test(Life_anova,trt = c("Treatment", "Life")) #HSD Tukey 
HSD.stat

## Coverage by Functional Life per Plot Box Plot #
Life_Box = 
  ggplot(LIFE, aes(x = Treatment, y = TotalCV, fill = Life)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(size = 2, alpha = 0.6, shape = 21, 
              position = position_jitterdodge(jitter.width = 0.1))
Life_Box

################################################################################
############################## Raunkiaer #######################################
################################################################################

RAUN <- GRIN %>%
  filter(Raunkier != "") %>%
  select(ID, Treatment, Raunkier, Species, Coverage) %>% 
  group_by(ID, Treatment, Raunkier)  %>% 
  summarize(TotalCV = sum(Coverage))

# Test for Significance #
Raunkier_anova = aov(TotalCV ~ Treatment + Raunkier, data = RAUN)
summary(Raunkier_anova)
tukey.one.way<-TukeyHSD(Raunkier_anova)
tukey.one.way
HSD.stat = HSD.test(Raunkier_anova,trt = c("Treatment", "Raunkier")) #HSD Tukey 
HSD.stat

## Coverage by Functional Raunkier per Plot Box Plot #
Raunkier_Box = 
  ggplot(RAUN, aes(x = Treatment, y = TotalCV, fill = Raunkier)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(size = 2, alpha = 0.6, shape = 21, 
              position = position_jitterdodge(jitter.width = 0.1))
Raunkier_Box

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

################################################################################
##################### Functional Groups Analysis & Plot ########################
################################################################################

GROUP <- GRIN %>% 
  filter(Group != "Bare") %>%
  select(ID, Treatment, Group, Species, Coverage) %>% 
  group_by(ID, Treatment, Group)  %>% 
  summarize(TotalCV = sum(Coverage))

# Test for Significance #
Fun_anova = aov(TotalCV ~ Treatment + Group, data = GROUP)
summary(Fun_anova)
tukey.one.way<-TukeyHSD(Fun_anova)
tukey.one.way
HSD.stat = HSD.test(Fun_anova,trt = c("Treatment", "Group")) #HSD Tukey 
HSD.stat

## Coverage by Functional Group per Plot Box Plot #
Group_Box = 
  ggplot(GROUP, aes(x = Treatment, y = TotalCV, fill = Group)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(size = 2, alpha = 0.6, shape = 21, 
              position = position_jitterdodge(jitter.width = 0.1))
Group_Box

################################################################################
############################## Annual vs Perennial #############################
################################################################################

LIFE <- GRIN %>% 
  filter(Life != "") %>%
  select(ID, Treatment, Life, Species, Coverage) %>% 
  group_by(ID, Treatment, Life)  %>% 
  summarize(TotalCV = sum(Coverage))

# Test for Significance #
Life_anova = aov(TotalCV ~ Treatment + Life, data = LIFE)
summary(Life_anova)
tukey.one.way<-TukeyHSD(Life_anova)
tukey.one.way
HSD.stat = HSD.test(Life_anova,trt = c("Treatment", "Life")) #HSD Tukey 
HSD.stat

## Coverage by Functional Life per Plot Box Plot #
Life_Box = 
  ggplot(LIFE, aes(x = Treatment, y = TotalCV, fill = Life)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(size = 2, alpha = 0.6, shape = 21, 
              position = position_jitterdodge(jitter.width = 0.1))
Life_Box

################################################################################
############################## Raunkiaer #######################################
################################################################################

RAUN <- GRIN %>%
  filter(Raunkier != "") %>%
  select(ID, Treatment, Raunkier, Species, Coverage) %>% 
  group_by(ID, Treatment, Raunkier)  %>% 
  summarize(TotalCV = sum(Coverage))

# Test for Significance #
Raunkier_anova = aov(TotalCV ~ Treatment + Raunkier, data = RAUN)
summary(Raunkier_anova)
tukey.one.way<-TukeyHSD(Raunkier_anova)
tukey.one.way
HSD.stat = HSD.test(Raunkier_anova,trt = c("Treatment", "Raunkier")) #HSD Tukey 
HSD.stat

## Coverage by Functional Raunkier per Plot Box Plot #
Raunkier_Box = 
  ggplot(RAUN, aes(x = Treatment, y = TotalCV, fill = Raunkier)) +
  geom_boxplot(outlier.shape = NA ) +
  geom_jitter(size = 2, alpha = 0.6, shape = 21, 
              position = position_jitterdodge(jitter.width = 0.1))
Raunkier_Box
