################################################################################
################################################################################
#########################    GRIN - Disturbance   ##############################
#########################   Species Coverage      ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2020 - 2022        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", 
                      "ggsignif", "multcompView", "ggpubr", "rstatix",
                      "vegan", "labdsv", "tables", "plotrix")
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
library(tables)
library(plotrix)

####################### Read in 2021 - 2023 Data  ##############################
Data = read.csv("Data/GRIN - 2020-2023.csv")
Data$Coverage = as.numeric(Data$Coverage)

str(Data)
summary(Data)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
Data <- mutate(Data, Coverage = case_when(
  grepl(10, Coverage) ~ 97.5,
  grepl(1, Coverage) ~ 0.1,
  grepl(2, Coverage) ~ 0.5,
  grepl(3, Coverage) ~ 1.5,
  grepl(4, Coverage) ~ 3.5,
  grepl(5, Coverage) ~ 7.5,
  grepl(6, Coverage) ~ 17.5,
  grepl(7, Coverage) ~ 37.5,
  grepl(8, Coverage) ~ 62.5,
  grepl(9, Coverage) ~ 85,
  grepl(0, Coverage) ~ 0,
))

# Remove Seeding Treatment # 
Data = filter(Data, Treatment != 'S')

################################################################################
################################################################################
################################################################################
############################### Hairy Indigo ###################################
################################################################################
################################################################################
################################################################################

Ih = filter(Data, Species == "Indigofera hirsuta")
summary(Ih)

# Creates data sets by year #
Ih_21 = filter(Ih, Year == 1)
Ih_22 = filter(Ih, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Ih_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ih_21$Treatment= as.factor(Ih_21$Treatment)
Ih_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_21 = Ih_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- Ih_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

############################### 2022 Ih ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Ih_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ih_22$Treatment= as.factor(Ih_22$Treatment)
Ih_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Ih_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Ih_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

## Lovegrass Coverage 2022 Box plot ##
IhBox21 = 
  ggplot(Ih_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "I. hirsuta % Coverage", title = "2021")
IhBox21

## Pityopsis trayci  Coverage 2023 Boxplot ##
IhBox22 = 
  ggplot(Ih_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "I. hirsuta % Coverage", title = "2022")
IhBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Ih_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ih_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Ih_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ih_22.csv")

################################################################################
################################################################################
################################################################################
############################### Oxalis #########################################
################################################################################
################################################################################
################################################################################

Oc = filter(Data, Species == "Oxalis corniculata")
summary(Oc)

# Creates data sets by year #
Oc_21 = filter(Oc, Year == 1)
Oc_22 = filter(Oc, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Oc_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Oc_21$Treatment= as.factor(Oc_21$Treatment)
Oc_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_21 = Oc_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- Oc_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

############################### 2022 Oc ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Oc_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Oc_22$Treatment= as.factor(Oc_22$Treatment)
Oc_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Oc_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Oc_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

## Lovegrass Coverage 2022 Box plot ##
OcBox21 = 
  ggplot(Oc_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "O. corniculata % Coverage", title = "2021")
OcBox21

## Pityopsis trayci  Coverage 2023 Boxplot ##
OcBox22 = 
  ggplot(Oc_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "O. corniculata % Coverage", title = "2022")
OcBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Oc_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Oc_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Oc_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Oc_22.csv")

################################################################################
################################################################################
################################################################################
############################### P. setaceum ####################################
################################################################################
################################################################################
################################################################################

Ps = filter(Data, Species == "Paspalum setaceum")
summary(Ps)

# Creates data sets by year #
Ps_21 = filter(Ps, Year == 1)
Ps_22 = filter(Ps, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Ps_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ps_21$Treatment= as.factor(Ps_21$Treatment)
Ps_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_21 = Ps_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- Ps_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

############################### 2022 Ps ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Ps_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ps_22$Treatment= as.factor(Ps_22$Treatment)
Ps_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Ps_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Ps_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

## Lovegrass Coverage 2022 Box plot ##
PsBox21 = 
  ggplot(Ps_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "P. setaceum % Coverage", title = "2021")
PsBox21

## Pityopsis trayci  Coverage 2023 Boxplot ##
PsBox22 = 
  ggplot(Ps_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "P. setaceum % Coverage", title = "2022")
PsBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Ps_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ps_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Ps_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ps_22.csv")


################################################################################
################################################################################
################################################################################
############################### Bahia ##########################################
################################################################################
################################################################################
################################################################################

Pn = filter(Data, Species == "Paspalum notatum")
summary(Pn)

# Creates data sets by year #
Pn_20 = filter(Pn, Year == 0)
Pn_21 = filter(Pn, Year == 1)
Pn_22 = filter(Pn, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2020 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Pn_20)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Pn_20$Treatment= as.factor(Pn_20$Treatment)
Pn_20 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_20 = Pn_20 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_20)

tukey_20 <- Pn_20 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_20

############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Pn_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Pn_21$Treatment= as.factor(Pn_21$Treatment)
Pn_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_21 = Pn_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- Pn_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

############################### 2022 Pn ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Pn_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Pn_22$Treatment= as.factor(Pn_22$Treatment)
Pn_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Pn_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Pn_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

## Coverage 2020 Boxplot ##
PnBox20 = 
  ggplot(Pn_20, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_20, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_20,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_20)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "P. notatum % Coverage", title = "2020")
PnBox20

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Pn_20)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_20.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Pn_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Pn_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_22.csv")

##  Coverage 2021 Box plot ##
PnBox21 = 
  ggplot(Pn_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "P. notatum % Coverage", title = "2021")
PnBox21

## Coverage 2023 Boxplot ##
PnBox22 = 
  ggplot(Pn_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "P. notatum % Coverage", title = "2022")
PnBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Pn_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Pn_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_22.csv")

################################################################################
################################################################################
################################################################################
############################### Richardia ######################################
################################################################################
################################################################################
################################################################################

Rs = filter(Data, Species == "Richardia spp.")
summary(Rs)

# Creates data sets by year #
Rs_21 = filter(Rs, Year == 1)
Rs_22 = filter(Rs, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Rs_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Rs_21$Treatment= as.factor(Rs_21$Treatment)
Rs_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_21 = Rs_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- Rs_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

############################### 2022 Rs ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Rs_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Rs_22$Treatment= as.factor(Rs_22$Treatment)
Rs_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Rs_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Rs_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

## Lovegrass Coverage 2022 Box plot ##
RsBox21 = 
  ggplot(Rs_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Richardia spp. % Coverage", title = "2021")
RsBox21

## Pityopsis trayci  Coverage 2023 Boxplot ##
RsBox22 = 
  ggplot(Rs_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Richardia spp. % Coverage", title = "2022")
RsBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Rs_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Rs_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Rs_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Rs_22.csv")

################################################################################
################################################################################
################################################################################
############################### Dichondra ######################################
################################################################################
################################################################################
################################################################################

Ds = filter(Data, Species == "Dichondra spp.")
summary(Ds)

# Creates data sets by year #
Ds_21 = filter(Ds, Year == 1)
Ds_22 = filter(Ds, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Ds_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ds_21$Treatment= as.factor(Ds_21$Treatment)
Ds_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_21 = Ds_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- Ds_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

############################### 2022 Ds ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Ds_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ds_22$Treatment= as.factor(Ds_22$Treatment)
Ds_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Ds_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Ds_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

##  Coverage 2022 Box plot ##
DsBox21 = 
  ggplot(Ds_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Dichondra spp. % Coverage", title = "2021")
DsBox21

## Coverage 2023 Boxplot ##
DsBox22 = 
  ggplot(Ds_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Dichondra spp. % Coverage", title = "2022")
DsBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Ds_21)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ds_21.csv")

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Ds_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ds_22.csv")


################################################################################
################################################################################
################################################################################
############################### Hairy Indigo ###################################
################################################################################
################################################################################
################################################################################

Eo = filter(Data, Species == "Eremochloa ophiuroides")
summary(Eo)

# Creates data sets by year #
Eo_21 = filter(Eo, Year == 1)
Eo_22 = filter(Eo, Year == 2)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Eo ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = Eo_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Eo_22$Treatment= as.factor(Eo_22$Treatment)
Eo_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_22 = Eo_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Eo_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

## Lovegrass Coverage 2022 Box plot ##
EoBox21 = 
  ggplot(Eo_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_21,size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "E. ophiuroides % Coverage", title = "2021")
EoBox21

## Coverage 2023 Boxplot ##
EoBox22 = 
  ggplot(Eo_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "E. ophiuroides % Coverage", title = "2022")
EoBox22

tmp <- tabular(Treatment ~ Coverage * (mean+sd+std.error), data=Eo_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Eo_22.csv")

################## Save Figures Above using ggarrange ##########################
Change = 
  ggarrange(IhBox21, OcBox21, PsBox21,PnBox21, RsBox21, DsBox21, 
            EoBox21, ncol = 2, nrow = 4)
Change
annotate_figure(Change, top = text_grob("", color = "black", 
                                        face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/species21.png", 
       width = 10, height = 18)

Change = 
  ggarrange(IhBox22, OcBox22, PsBox22, PnBox22, RsBox22, DsBox22, 
            EoBox22, ncol = 2, nrow = 4)
Change
annotate_figure(Change, top = text_grob("", color = "black", 
                                        face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/species22.png", 
       width = 10, height = 18)
