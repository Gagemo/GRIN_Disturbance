################################################################################
################################################################################
######################           GRIN - Fire                ####################
######################        Species Richness              ####################
######################      University of Florida           ####################
######################          Gage LaPierre               ####################
######################           2020 - 2022                ####################
################################################################################
################################################################################

######################### Clears Environment & History  ########################

rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################

list.of.packages <- c( "showtext", "tidyverse", "vegan", "labdsv", "reshape", 
                       "ggpubr", "rstatix")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(reshape)
library(showtext)
library(ggpubr)
library(rstatix)

##########################     Read in Data       ##############################
data = read.csv("Data/GRIN - 2020-2023.csv")
data$Coverage = as.numeric(data$Coverage)

str(data)
summary(data)

# Select just Seeding Treatment # 
data = filter(data, Treatment == 'S')

### Remove NA and empty Values ###
data = filter(data, Coverage != "NA") %>%
  filter(Coverage != "")

# Creates data sets by year #
data_22 = filter(data, Year == 2)
data_23 = filter(data, Year == 3)

# Create Species Pivot Table by Year #
Spp_22 = dplyr::select(data_22, ID, Species, Coverage) %>% matrify() 
Spp_22[] <- lapply(Spp_22, as.numeric)

Spp_23 = dplyr::select(data_23, ID, Species, Coverage) %>% matrify() 
Spp_23[] <- lapply(Spp_23, as.numeric)

# Calculates species richness by year #
SR_22 = specnumber(Spp_22)
SR_22 = as.data.frame(SR_22)

SR_23 = specnumber(Spp_23)
SR_23 = as.data.frame(SR_23)

# Create Grouped Treatment Table and Summaries to fit Species Table #
Treat_22 <- select(data_22, Fire, ID) %>% 
  group_by(ID, Fire) %>% 
  dplyr::summarise()

Treat_23 <- select(data_23, Fire, ID) %>% 
  group_by(ID, Fire) %>% 
  dplyr::summarise()

## Merge species richness with habitat/plot data for ggplot ##
SR_treat_22 = cbind(Treat_22, SR_22)
colnames(SR_treat_22)[3] <- "SR"
as.numeric(SR_treat_22$SR)
SR_treat_22$SR = as.numeric(SR_treat_22$SR)
SR_treat_22 = as.data.frame(SR_treat_22)

SR_treat_23 = cbind(Treat_23, SR_23)
colnames(SR_treat_23)[3] <- "SR"
SR_treat_23$SR = as.numeric(SR_treat_23$SR)
SR_treat_23 = as.data.frame(SR_treat_23)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(SR ~ Fire, data = SR_treat_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
SR_treat_22$Fire = as.factor(SR_treat_22$Fire)
SR_treat_22 %>% levene_test(SR_treat_22$SR ~ SR_treat_22$Fire)

# Test for Significance #
anova_22 = SR_treat_22 %>% anova_test(SR ~ Fire) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- SR_treat_22 %>% 
  tukey_hsd(SR ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Data ######################################
# Check Assumptions #
model  <- lm(SR ~ Fire, data = SR_treat_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
SR_treat_23$Fire = as.factor(SR_treat_23$Fire)
SR_treat_23 %>% levene_test(SR_treat_23$SR ~ SR_treat_23$Fire)
anova_23 = SR_treat_23 %>% anova_test(SR ~ Fire) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- SR_treat_23 %>% 
  tukey_hsd(SR ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

################################################################################
################ Create SR Boxplots for each year ##############################
################################################################################

SR_Box_22 = 
  ggplot(SR_treat_22, aes(x = Fire, y = SR), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), size = 4, alpha = 0.7) +
  stat_pvalue_manual(tukey_22,
                     hide.ns = T)+
  ylim(0, 30) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
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
        text=element_text(size=22),
        axis.title.x = element_text(size=22, face="bold", colour = "black"),    
        axis.title.y = element_text(size=22, face="italic", colour = "black"),   
        axis.text.x=element_text(size=22, face = "bold", color = "black"),
        axis.text.y=element_text(size=22, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 22, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Species Richness", title = "2022")
SR_Box_22
ggsave("Figures/Chapter 2 - Fire/SR_Box_2022.png")

SR_Box_23 = 
  ggplot(SR_treat_23, aes(x = Fire, y = SR), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), size = 5, alpha = 0.7) +
  stat_pvalue_manual(tukey_23,
                     hide.ns = T)+
  ylim(0, 30) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
        text=element_text(size=22),
        axis.title.x = element_text(size=22, face="bold", colour = "black"),    
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=22, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 22, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Species Richness", title = "2023")
SR_Box_23
ggsave("Figures/Chapter 2 - Fire/SR_Box_2023.png")

################## Save Figures Above using ggarrange ##########################
ggarrange(SR_Box_22, SR_Box_23, ncol = 2, nrow = 1)
ggsave("Figures/Chapter 2 - Fire/2022-2023_SR.png", 
       width = 16, height = 10)

