################################################################################
################################################################################
#########################       GRIN - Fire       ##############################
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
                      "vegan", "labdsv")
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

##########################     Read in 2020 - 2023 Data       ##################
GRIN = read.csv("Data/GRIN_FUN-2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Select just Seeding Treatment # 
GRIN = filter(GRIN, Treatment == 'S')

# Remove Bareground Treatment # 
GRIN = filter(GRIN, Group != 'Bare')

# Replace NA Values with Zeros#
GRIN$Coverage[is.na(GRIN$Coverage)] <- 0

GRIN$Coverage = as.numeric(GRIN$Coverage)

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

# Creates data sets by year #
GRIN_22 = filter(GRIN, Year == 2)
GRIN_23 = filter(GRIN, Year == 3)

################################################################################
####################### Functional Groups Analysis #############################
################################################################################
########################### 2022 Data Set ######################################
################################################################################

# Sums CV data #
GRIN_22M <- dplyr::select(GRIN_22, ID, Group, Coverage) %>% matrify()
GRIN_22M <- tibble::rownames_to_column(GRIN_22M, "ID")

# Create grouped treatment/environment table and summaries to fit species table#
Treat_22 = group_by(GRIN_22, ID, Fire) %>% dplyr::summarize()

###############################  Forbs 2022 ####################################
Forb = select(GRIN_22M, Forb) %>%
  cbind(Treat_22)
Forb <- ungroup(Forb)

# Check Assumptions #
model  <- lm(Forb ~ Fire, data = Forb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb$Fire = as.factor(Forb$Fire)
Forb %>% levene_test(Forb ~ Fire)

# Test for Significance #
anova_Forb = Forb %>% anova_test(Forb ~ Fire) %>% 
  add_significance()
summary(anova_Forb)

tukey_Forb <- Forb %>% 
  tukey_hsd(Forb ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Forb

###############################  Grass 2022 ####################################
Grass = select(GRIN_22M, Grass) %>%
  cbind(Treat_22)
Grass <- ungroup(Grass)

# Check Assumptions #
model  <- lm(Grass ~ Fire, data = Grass)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Grass$Fire = as.factor(Grass$Fire)
Grass %>% levene_test(Grass ~ Fire)

# Test for Significance #
anova_Grass = Grass %>% anova_test(Grass ~ Fire) %>% 
  add_significance()
summary(anova_Grass)

tukey_Grass <- Grass %>% 
  tukey_hsd(Grass ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Grass


###############################  Sedges 2022 ###################################
Sedge = select(GRIN_22M, Sedge) %>%
  cbind(Treat_22)
Sedge <- ungroup(Sedge)

# Check Assumptions #
model  <- lm(Sedge ~ Fire, data = Sedge)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Sedge$Fire = as.factor(Sedge$Fire)
Sedge %>% levene_test(Sedge ~ Fire)

# Test for Significance #
anova_Sedge = Sedge %>% kruskal_test(Sedge ~ Fire) %>% 
  add_significance()
summary(anova_Sedge)

tukey_Sedge <- Sedge %>% 
  dunn_test(Sedge ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Sedge

###############################  Woody 2022  ###################################
Woody = select(GRIN_22M, Woody) %>%
  cbind(Treat_22)
Woody <- ungroup(Woody)

# Check Assumptions #
model  <- lm(Woody ~ Fire, data = Woody)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Woody$Fire = as.factor(Woody$Fire)
Woody %>% levene_test(Woody ~ Fire)

# Test for Significance #
anova_Woody = Woody %>% kruskal_test(Woody ~ Fire) %>% 
  add_significance()
summary(anova_Woody)

tukey_Woody <- Woody %>% 
  dunn_test(Woody ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Woody

################################################################################
########### Coverage by Functional Group per Plot Box Plot #####################
################################################################################

# Forbs 2022 #
Forb_Box = 
  ggplot(Forb, aes(x = Fire, y = Forb), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Forb,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Forb,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Forb)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Forb")
Forb_Box

# Grasses 2022 #
Grass_Box = 
  ggplot(Grass, aes(x = Fire, y = Grass), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Grass,size = 8, bracket.size = 1, hide.ns = T) +
  labs(subtitle = get_test_label(anova_Grass,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Grass)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Grass")
Grass_Box

# Sedges 2022 #
Sedge_Box = 
  ggplot(Sedge, aes(x = Fire, y = Sedge), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Sedge,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Sedge,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Sedge)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Sedge")
Sedge_Box

# Woody 2022 #
Woody_Box = 
  ggplot(Woody, aes(x = Fire, y = Woody), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Woody,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Woody,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Woody)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Woody")
Woody_Box

################## Save Figures Above using ggarrange ##########################
CombinedFun_2022 = 
  ggarrange(Forb_Box, Grass_Box, Sedge_Box, Woody_Box, ncol = 2, nrow = 2)
annotate_figure(CombinedFun_2022, top = text_grob("2022", color = "black", 
                                                  face = "bold", size = 25))
ggsave("Figures/Chapter 2 - Fire/2022_Fun.png", 
       width = 12, height = 8)

################################################################################
####################### Functional Groups Analysis #############################
################################################################################
########################### 2023 Data Set ######################################
################################################################################

# Sums CV data #
GRIN_23M <- dplyr::select(GRIN_23, ID, Group, Coverage) %>% matrify()
GRIN_23M <- tibble::rownames_to_column(GRIN_23M, "ID")

# Create grouped treatment/environment table and summaries to fit species table#
Treat_23 = group_by(GRIN_23, ID, Fire) %>% dplyr::summarize()

###############################  Forbs 2023 ####################################
Forb = select(GRIN_23M, Forb) %>%
  cbind(Treat_23)
Forb <- ungroup(Forb)

# Check Assumptions #
model  <- lm(Forb ~ Fire, data = Forb)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb$Fire = as.factor(Forb$Fire)
Forb %>% levene_test(Forb ~ Fire)

# Test for Significance #
anova_Forb = Forb %>% kruskal_test(Forb ~ Fire) %>% 
  add_significance()
summary(anova_Forb)

tukey_Forb <- Forb %>% 
  dunn_test(Forb ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Forb

###############################  Grass 2023 ####################################
Grass = select(GRIN_23M, Grass) %>%
  cbind(Treat_23)
Grass <- ungroup(Grass)

# Check Assumptions #
model  <- lm(Grass ~ Fire, data = Grass)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Grass$Fire = as.factor(Grass$Fire)
Grass %>% levene_test(Grass ~ Fire)

# Test for Significance #
anova_Grass = Grass %>% anova_test(Grass ~ Fire) %>% 
  add_significance()
summary(anova_Grass)

tukey_Grass <- Grass %>% 
  tukey_hsd(Grass ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Grass

###############################  Sedges 2023 ###################################
Sedge = select(GRIN_23M, Sedge) %>%
  cbind(Treat_23)
Sedge <- ungroup(Sedge)

# Check Assumptions #
model  <- lm(Sedge ~ Fire, data = Sedge)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Sedge$Fire = as.factor(Sedge$Fire)
Sedge %>% levene_test(Sedge ~ Fire)

# Test for Significance #
anova_Sedge = Sedge %>% kruskal_test(Sedge ~ Fire) %>% 
  add_significance()
summary(anova_Sedge)

tukey_Sedge <- Sedge %>% 
  dunn_test(Sedge ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Sedge

###############################  Woody 2023  ###################################
Woody = select(GRIN_23M, Woody) %>%
  cbind(Treat_23)
Woody <- ungroup(Woody)

# Check Assumptions #
model  <- lm(Woody ~ Fire, data = Woody)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Woody$Fire = as.factor(Woody$Fire)
Woody %>% levene_test(Woody ~ Fire)

# Test for Significance #
anova_Woody = Woody %>% anova_test(Woody ~ Fire) %>% 
  add_significance()
summary(anova_Woody)

tukey_Woody <- Woody %>% 
  tukey_hsd(Woody ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Woody

################################################################################
########### Coverage by Functional Group per Plot Box Plot #####################
################################################################################

# Forbs 2023 #
Forb_Box = 
  ggplot(Forb, aes(x = Fire, y = Forb), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Forb,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Forb,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Forb)) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Forb")
Forb_Box
ggsave("Figures/Chapter 2 - Fire/2023Forb.png", 
       width = 10, height = 7)

# Grasses 2023 #
Grass_Box = 
  ggplot(Grass, aes(x = Fire, y = Grass), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Grass,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Grass,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Grass)) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Grass")
Grass_Box
ggsave("Figures/Chapter 2 - Fire/2023Grass_Box.png", 
       width = 10, height = 7)

# Sedges 2023 #
Sedge_Box = 
  ggplot(Sedge, aes(x = Fire, y = Sedge), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Sedge,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Sedge,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Sedge)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
  ylim(0, 17) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Sedge")
Sedge_Box
ggsave("Figures/Chapter 2 - Fire/2023Sedge_Box.png", 
       width = 10, height = 7)

# Woody 2023 #
Woody_Box = 
  ggplot(Woody, aes(x = Fire, y = Woody), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), alpha = 0.5) +
  stat_pvalue_manual(tukey_Woody,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Woody,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_Woody)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No Burn', 'Late-Spring', 'Winter')) +
  ylim(0,90) +
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
  labs(x = "", y = "Total (%) Coverage", title = "Woody")
Woody_Box
ggsave("Figures/Chapter 2 - Fire/2023Woody_Box.png", 
       width = 10, height = 7)

################## Save Figures Above using ggarrange ##########################
CombinedFun_2023 = 
  ggarrange(Forb_Box, Grass_Box, Sedge_Box, Woody_Box, ncol = 2, nrow = 2)
annotate_figure(CombinedFun_2023, top = text_grob("2023", color = "black", 
                                                  face = "bold", size = 25))
ggsave("Figures/Chapter 2 - Fire/2023_Fun.png", 
       width = 12, height = 8)

