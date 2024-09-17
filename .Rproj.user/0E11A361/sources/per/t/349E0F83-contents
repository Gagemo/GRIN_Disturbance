################################################################################
################################################################################
#########################   GRIN - Disturbance    ##############################
#########################    Functional Groups    ##############################
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
                      "vegan", "labdsv", "tables")
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

##########################     Read in 2020 - 2023 Data       ##################
GRIN = read.csv("Data/GRIN - 2020-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Year != 3)

# Remove Bareground Treatment # 
GRIN = filter(GRIN, Group != 'Bare')

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

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
  grepl(9, Coverage) ~ 85,
  
))

# Creates data sets by year #
GRIN_20 = filter(GRIN, Year == 0)
GRIN_21 = filter(GRIN, Year == 1)
GRIN_22 = filter(GRIN, Year == 2)


################################################################################
####################### Functional Groups Analysis #############################
################################################################################
########################### 2020 Data Set ######################################
################################################################################

# Sums CV data #
GRIN_20M <- dplyr::select(GRIN_20, ID, Group, Coverage) %>% matrify()
GRIN_20M <- tibble::rownames_to_column(GRIN_20M, "ID")

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_20, ID, Treatment) %>% dplyr::summarize()

###############################  Forbs 2020 ####################################
Forb = select(GRIN_20M, Forb) %>%
  cbind(Treat)
Forb <- ungroup(Forb)

# Check Assumptions #
model  <- lm(Forb ~ Treatment, data = Forb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb$Treatment = as.factor(Forb$Treatment)
Forb %>% levene_test(Forb ~ Treatment)

# Test for Significance #
anova_Forb = Forb %>% kruskal_test(Forb ~ Treatment) %>% 
  add_significance()
summary(anova_Forb)

tukey_Forb <- Forb %>% 
  dunn_test(Forb ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Forb

tmp <- tabular(Treatment ~ Forb* (mean+sd+std.error), data=Forb)
tmp

###############################  Grass 2020 ####################################
Grass = select(GRIN_20M, Grass) %>%
  cbind(Treat)
Grass <- ungroup(Grass)

# Check Assumptions #
model  <- lm(Grass ~ Treatment, data = Grass)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Grass$Treatment = as.factor(Grass$Treatment)
Grass %>% levene_test(Grass ~ Treatment)

# Test for Significance #
anova_Grass = Grass %>% kruskal_test(Grass ~ Treatment) %>% 
  add_significance()
summary(anova_Grass)

tukey_Grass <- Grass %>% 
  dunn_test(Grass ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Grass
tmp <- tabular(Treatment ~ Grass* (mean+sd+std.error), data=Grass)
tmp

###############################  Sedges 2020 ###################################
Sedge = select(GRIN_20M, Sedge) %>%
  cbind(Treat)
Sedge <- ungroup(Sedge)

# Check Assumptions #
model  <- lm(Sedge ~ Treatment, data = Sedge)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Sedge$Treatment = as.factor(Sedge$Treatment)
Sedge %>% levene_test(Sedge ~ Treatment)

# Test for Significance #
anova_Sedge = Sedge %>% kruskal_test(Sedge ~ Treatment) %>% 
  add_significance()
summary(anova_Sedge)

tukey_Sedge <- Sedge %>% 
  dunn_test(Sedge ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Sedge

tmp <- tabular(Treatment ~ Sedge* (mean+sd+std.error), data=Sedge)
tmp

################################################################################
########### Coverage by Functional Group per Plot Box Plot #####################
################################################################################

# Forbs 2020 #
Forb_Box = 
  ggplot(Forb, aes(x = Treatment, y = Forb), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Forb, size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Forb,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Forb)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Forb")
Forb_Box

# Grasses 2020 #
Grass_Box = 
  ggplot(Grass, aes(x = Treatment, y = Grass), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Grass,hide.ns = T, size = 8, bracket.size = 1) +
  labs(subtitle = get_test_label(anova_Grass,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Grass)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Grass")
Grass_Box

# Sedges 2020 #
Sedge_Box = 
  ggplot(Sedge, aes(x = Treatment, y = Sedge), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Sedge,size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Sedge,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Sedge)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Sedge")
Sedge_Box

################## Save Figures Above using ggarrange ##########################
CombinedFun_2020 = 
  ggarrange(Forb_Box, Grass_Box, Sedge_Box, ncol = 2, nrow = 2)
annotate_figure(CombinedFun_2020, top = text_grob("2020", color = "black", 
                                                  face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2020_Fun.png", 
       width = 12, height = 8)


################################################################################
####################### Functional Groups Analysis #############################
################################################################################
########################### 2021 Data Set ######################################
################################################################################

# Sums CV data #
GRIN_21M <- dplyr::select(GRIN_21, ID, Group, Coverage) %>% matrify()
GRIN_21M <- tibble::rownames_to_column(GRIN_21M, "ID")

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_21, ID, Treatment) %>% dplyr::summarize()

###############################  Forbs 2021 ####################################
Forb = select(GRIN_21M, Forb) %>%
  cbind(Treat)
Forb <- ungroup(Forb)

# Check Assumptions #
model  <- lm(Forb ~ Treatment, data = Forb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb$Treatment = as.factor(Forb$Treatment)
Forb %>% levene_test(Forb ~ Treatment)

# Test for Significance #
anova_Forb = Forb %>% kruskal_test(Forb ~ Treatment) %>% 
  add_significance()
summary(anova_Forb)

tukey_Forb <- Forb %>% 
  dunn_test(Forb ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Forb

tmp <- tabular(Treatment ~ Forb* (mean+sd+std.error), data=Forb)
tmp

###############################  Grass 2021 ####################################
Grass = select(GRIN_21M, Grass) %>%
  cbind(Treat)
Grass <- ungroup(Grass)

# Check Assumptions #
model  <- lm(Grass ~ Treatment, data = Grass)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Grass$Treatment = as.factor(Grass$Treatment)
Grass %>% levene_test(Grass ~ Treatment)

# Test for Significance #
anova_Grass = Grass %>% anova_test(Grass ~ Treatment) %>% 
  add_significance()
summary(anova_Grass)

tukey_Grass <- Grass %>% 
  tukey_hsd(Grass ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Grass
tmp <- tabular(Treatment ~ Grass* (mean+sd+std.error), data=Grass)
tmp

###############################  Sedges 2021 ###################################
Sedge = select(GRIN_21M, Sedge) %>%
  cbind(Treat)
Sedge <- ungroup(Sedge)

# Check Assumptions #
model  <- lm(Sedge ~ Treatment, data = Sedge)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Sedge$Treatment = as.factor(Sedge$Treatment)
Sedge %>% levene_test(Sedge ~ Treatment)

# Test for Significance #
anova_Sedge = Sedge %>% kruskal_test(Sedge ~ Treatment) %>% 
  add_significance()
summary(anova_Sedge)

tukey_Sedge <- Sedge %>% 
  dunn_test(Sedge ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Sedge

tmp <- tabular(Treatment ~ Sedge* (mean+sd+std.error), data=Sedge)
tmp

###############################  Woody 2021  ###################################
Woody = select(GRIN_21M, Woody) %>%
  cbind(Treat)
Woody <- ungroup(Woody)

# Check Assumptions #
model  <- lm(Woody ~ Treatment, data = Woody)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Woody$Treatment = as.factor(Woody$Treatment)
Woody %>% levene_test(Woody ~ Treatment)

# Test for Significance #
anova_Woody = Woody %>% kruskal_test(Woody ~ Treatment) %>% 
  add_significance()
summary(anova_Woody)

tukey_Woody <- Woody %>% 
  dunn_test(Woody ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Woody

tmp <- tabular(Treatment ~ Woody* (mean+sd+std.error), data=Woody)
tmp

################################################################################
########### Coverage by Functional Group per Plot Box Plot #####################
################################################################################

# Forbs 2021 #
Forb_Box = 
  ggplot(Forb, aes(x = Treatment, y = Forb), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Forb, size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Forb,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Forb)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Forb")
Forb_Box

# Grasses 2021 #
Grass_Box = 
  ggplot(Grass, aes(x = Treatment, y = Grass), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Grass,hide.ns = T, size = 8, bracket.size = 1) +
  labs(subtitle = get_test_label(anova_Grass,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Grass)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Grass")
Grass_Box

# Sedges 2021 #
Sedge_Box = 
  ggplot(Sedge, aes(x = Treatment, y = Sedge), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Sedge,size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Sedge,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Sedge)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Sedge")
Sedge_Box

# Woody 2021 #
Woody_Box = 
  ggplot(Woody, aes(x = Treatment, y = Woody), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Woody, size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Woody,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Woody)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Woody")
Woody_Box

################## Save Figures Above using ggarrange ##########################
CombinedFun_2021 = 
  ggarrange(Forb_Box, Grass_Box, Sedge_Box, Woody_Box, ncol = 2, nrow = 2)
annotate_figure(CombinedFun_2021, top = text_grob("2021", color = "black", 
                                                  face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2021_Fun.png", 
       width = 12, height = 8)

################################################################################
####################### Functional Groups Analysis #############################
################################################################################
########################### 2022 Data Set ######################################
################################################################################

# Sums CV data #
GRIN_22M <- dplyr::select(GRIN_22, ID, Group, Coverage) %>% matrify()
GRIN_22M <- tibble::rownames_to_column(GRIN_22M, "ID")

# Create grouped treatment/environment table and summaries to fit species table#
Treat = group_by(GRIN_22, ID, Treatment) %>% dplyr::summarize()

###############################  Forbs 2022 ####################################
Forb = select(GRIN_22M, Forb) %>%
  cbind(Treat)
Forb <- ungroup(Forb)

# Check Assumptions #
model  <- lm(Forb ~ Treatment, data = Forb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb$Treatment = as.factor(Forb$Treatment)
Forb %>% levene_test(Forb ~ Treatment)

# Test for Significance #
anova_Forb = Forb %>% kruskal_test(Forb ~ Treatment) %>% 
  add_significance()
summary(anova_Forb)

tukey_Forb <- Forb %>% 
  dunn_test(Forb ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Forb

tmp <- tabular(Treatment ~ Forb* (mean+sd+std.error), data=Forb)
tmp

###############################  Grass 2022 ####################################
Grass = select(GRIN_22M, Grass) %>%
  cbind(Treat)
Grass <- ungroup(Grass)

# Check Assumptions #
model  <- lm(Grass ~ Treatment, data = Grass)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Grass$Treatment = as.factor(Grass$Treatment)
Grass %>% levene_test(Grass ~ Treatment)

# Test for Significance #
anova_Grass = Grass %>% anova_test(Grass ~ Treatment) %>% 
  add_significance()
summary(anova_Grass)

tukey_Grass <- Grass %>% 
  tukey_hsd(Grass ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Grass

tmp <- tabular(Treatment ~ Grass* (mean+sd+std.error), data=Grass)
tmp

###############################  Sedges 2022 ###################################
Sedge = select(GRIN_22M, Sedge) %>%
  cbind(Treat)
Sedge <- ungroup(Sedge)

# Check Assumptions #
model  <- lm(Sedge ~ Treatment, data = Sedge)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Sedge$Treatment = as.factor(Sedge$Treatment)
Sedge %>% levene_test(Sedge ~ Treatment)

# Test for Significance #
anova_Sedge = Sedge %>% kruskal_test(Sedge ~ Treatment) %>% 
  add_significance()
summary(anova_Sedge)

tukey_Sedge <- Sedge %>% 
  dunn_test(Sedge ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Sedge

tmp <- tabular(Treatment ~ Sedge* (mean+sd+std.error), data=Sedge)
tmp

###############################  Woody 2022  ###################################
Woody = select(GRIN_22M, Woody) %>%
  cbind(Treat)
Woody <- ungroup(Woody)

# Check Assumptions #
model  <- lm(Woody ~ Treatment, data = Woody)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Woody$Treatment = as.factor(Woody$Treatment)
Woody %>% levene_test(Woody ~ Treatment)

# Test for Significance #
anova_Woody = Woody %>% kruskal_test(Woody ~ Treatment) %>% 
  add_significance()
summary(anova_Woody)

tukey_Woody <- Woody %>% 
  dunn_test(Woody ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Woody

tmp <- tabular(Treatment ~ Woody* (mean+sd+std.error), data=Woody)
tmp

################################################################################
########### Coverage by Functional Group per Plot Box Plot #####################
################################################################################

# Forbs 2022 #
Forb_Box = 
  ggplot(Forb, aes(x = Treatment, y = Forb), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Forb,size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Forb,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Forb)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Forb")
Forb_Box

# Grasses 2022 #
Grass_Box = 
  ggplot(Grass, aes(x = Treatment, y = Grass), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Grass,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Grass,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Grass)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Grass")
Grass_Box

# Sedges 2022 #
Sedge_Box = 
  ggplot(Sedge, aes(x = Treatment, y = Sedge), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Sedge,size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Sedge,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Sedge)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Sedge")
Sedge_Box

# Woody 2022 #
Woody_Box = 
  ggplot(Woody, aes(x = Treatment, y = Woody), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Woody,size = 8, bracket.size = 1,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_Woody,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Woody)) +
  scale_fill_manual(labels=c('No-Till', 'Late-Spring', 'Winter'),
                    values=c("#FF3399", "#FFFF00", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
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
  labs(x = "Treatment", y = "Total (%) Coverage", title = "Woody")
Woody_Box

################## Save Figures Above using ggarrange ##########################
CombinedFun_2022 = 
  ggarrange(Forb_Box, Grass_Box, Sedge_Box, Woody_Box, ncol = 2, nrow = 2)
annotate_figure(CombinedFun_2022, top = text_grob("2022", color = "black", 
                                                  face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2022_Fun.png", 
       width = 12, height = 8)
