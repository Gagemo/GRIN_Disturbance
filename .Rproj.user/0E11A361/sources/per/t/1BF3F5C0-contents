################################################################################
################################################################################
#########################   Data - Disturbance    ##############################
#########################    Change in Cover      ##############################
#########################    Functional Groups    ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "extrafont", "plotrix", 
                      "ggsignif", "multcompView", "ggpubr", "rstatix", "labdsv",
                      "tables")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(agricolae)
library(extrafont)
library(ggsignif)
library(multcompView)
library(ggpubr)
library(plotrix)
library(rstatix)
library(tables)

##########################     Read in  Data       #############################
Data = read.csv("Data/GRIN - 2020-2023.csv")

Data$Coverage = as.numeric(Data$Coverage)

# Reclasifys coverage data (CV) from 1-10 scale to percent scale #
Data <- mutate(Data, Coverage = case_when(
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

str(Data)
summary(Data)

# Remove Seeding Treatment # 
Data = filter(Data, Treatment != 'S')

Data$YID <- paste(Data$Year,Data$ID)
Data$ID_ <- paste(Data$Treatment, Data$ID)

# Orders years and treatments so that they display in same sequence in graphs #
Data$Year = factor(Data$Year, levels=c('1','2'))

#################### Species abundances ########################################
# Creates and joins  data year 22 & 23 to make long data format #
Two_Abundance <- Data[which(Data$Year == "1"),]
Three_Abundance <- Data[which(Data$Year == "2"),]

Abundance_w <- full_join(Two_Abundance, Three_Abundance, 
                         by = c('ID_', "Treatment", 'Group','Species'))
Abundance_w = arrange(Abundance_w, Treatment)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(ID_, Treatment, Species, Group, 
                Coverage.x, Coverage.y) %>%
  group_by(ID_, Treatment, Species, Group) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x)

Change_Abundance = Change_Abundance %>%
  group_by(ID_, Treatment, Group) %>%
  summarise(total = sum(Change_abundance))

############################# WOODY COVER CAHNGES ##############################
woody = 
  Change_Abundance[which(Change_Abundance$Group == "Woody"),]
woody<-as.data.frame(woody)
woody$Treatment<-factor(woody$Treatment)

# Check Assumptions #
model  <- lm(total ~ Treatment, data = woody)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_woody = woody %>% anova_test(total ~ Treatment) %>% 
  add_significance()
anova_woody

lm(formula = total ~ Treatment, woody)
tukey_woody <- woody %>% 
  tukey_hsd(total ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_woody

tmp <- tabular(Treatment ~ total* (mean+sd+std.error), data=woody)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Woody_Change.csv")

woody_change_Box = 
  ggplot(woody, aes(x = Treatment, y = total), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_woody,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_woody, detailed = TRUE),
       caption = get_pwc_label(tukey_woody)) +
  ylim(-30,80) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Change in Coverage", title = "Woody Ruderals")
woody_change_Box

############################## FORB COVER CAHNGES ##############################
forb = 
  Change_Abundance[which(Change_Abundance$Group == "Forb"),]
forb<-as.data.frame(forb)
forb$Treatment<-factor(forb$Treatment)

# Check Assumptions #
model  <- lm(total ~ Treatment, data = forb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_forb = forb %>% anova_test(total ~ Treatment) %>% 
  add_significance()
anova_forb

lm(formula = total ~ Treatment, forb)
tukey_forb <- forb %>% 
  tukey_hsd(total ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_forb

tmp <- tabular(Treatment ~ total* (mean+sd+std.error), data=forb)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/forb_Change.csv")

forb_change_Box = 
  ggplot(forb, aes(x = Treatment, y = total), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_forb,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_forb, detailed = TRUE),
       caption = get_pwc_label(tukey_forb)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Change in Coverage", title = "Forbs")
forb_change_Box

#################################  Grass COVER CAHNGES ##############################
grass = 
  Change_Abundance[which(Change_Abundance$Group == "Grass"),]
grass<-as.data.frame(grass)
grass$Treatment<-factor(grass$Treatment)

# Check Assumptions #
model  <- lm(total ~ Treatment, data = grass)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_grass = grass %>% anova_test(total ~ Treatment) %>% 
  add_significance()
anova_grass

lm(formula = total ~ Treatment, grass)
tukey_grass <- grass %>% 
  tukey_hsd(total ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_grass

tmp <- tabular(Treatment ~ total* (mean+sd+std.error), data=grass)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/grass_Change.csv")

grass_change_Box = 
  ggplot(grass, aes(x = Treatment, y = total), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_grass,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_grass, detailed = TRUE),
       caption = get_pwc_label(tukey_grass)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Change in Coverage", title = "Grasses")
grass_change_Box

#################################  Sedge COVER CAHNGES ##############################
sedge = 
  Change_Abundance[which(Change_Abundance$Group == "Sedge"),]
sedge<-as.data.frame(sedge)
sedge$Treatment<-factor(sedge$Treatment)

# Check Assumptions #
model  <- lm(total ~ Treatment, data = sedge)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_sedge = sedge %>% anova_test(total ~ Treatment) %>% 
  add_significance()
anova_sedge

lm(formula = total ~ Treatment, sedge)
tukey_sedge <- sedge %>% 
  tukey_hsd(total ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_sedge

tmp <- tabular(Treatment ~ total* (mean+sd+std.error), data=sedge)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/sedge_Change.csv")

sedge_change_Box = 
  ggplot(sedge, aes(x = Treatment, y = total), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_sedge,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_sedge, detailed = TRUE),
       caption = get_pwc_label(tukey_sedge)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Change in Coverage", title = "Sedges")
sedge_change_Box

################## Save Figures Above using ggarrange ##########################
Change = 
  ggarrange(grass_change_Box, forb_change_Box, 
            sedge_change_Box, woody_change_Box, ncol = 2, nrow = 2)
annotate_figure(Change, top = text_grob("", color = "black", 
                                        face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/Changefun.png", 
       width = 12, height = 8)

######################## BAREGROUND COVER CAHNGES ##############################

bare = 
  Change_Abundance[which(Change_Abundance$Group == "Bare"),]
bare<-as.data.frame(bare)
bare$Treatment<-factor(bare$Treatment)

# Check Assumptions #
model  <- lm(total ~ Treatment, data = bare)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_bare = bare %>% kruskal_test(total ~ Treatment) %>% 
  add_significance()
anova_bare

lm(formula = total ~ Treatment, bare)
tukey_bare <- bare %>% 
  dunn_test(total ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_bare

tmp <- tabular(Treatment ~ total * (mean+sd+std.error), data=bare)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Bare_Change.csv")

bare_change_Box = 
  ggplot(bare, aes(x = Treatment, y = total), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 4, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_bare,size = 12, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_bare, detailed = TRUE),
       caption = get_pwc_label(tukey_bare)) +
  ylim(-10,100) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=20),
        axis.title.x = element_text(size=20, face="bold", colour = "black"),    
        axis.title.y = element_text(size=20, face="bold", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Change in Coverage", title = "Bare Ground")
bare_change_Box
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/change_Bareground.png", 
       width = 10, height = 7)
