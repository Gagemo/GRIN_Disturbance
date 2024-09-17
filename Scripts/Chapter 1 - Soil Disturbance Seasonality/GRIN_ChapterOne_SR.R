################################################################################
################################################################################
######################  GRIN - Soil Disturbance Seasonality ####################
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

list.of.packages <- c( "showtext", "tidyverse", "vegan", "labdsv", "reshape")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(labdsv)
library(reshape)
library(showtext)

##########################     Read in Data       ##############################
data = read.csv("Data/GRIN - 2020-2023.csv")
data$Coverage = as.numeric(data$Coverage)

str(data)
summary(data)

# Remove Seeding Treatment # 
data = filter(data, Treatment != 'S')

### Remove NA and empty Values ###
data = filter(data, Coverage != "NA") %>%
  filter(Coverage != "")

# Creates data sets by year #
data_20 = filter(data, Year == 0)
data_21 = filter(data, Year == 1)
data_22 = filter(data, Year == 2)

# Create Species Pivot Table by Year #
Spp_20 = dplyr::select(data_20, ID, Species, Coverage) %>% matrify() 
Spp_20[] <- lapply(Spp_20, as.numeric)

Spp_21 = dplyr::select(data_21, ID, Species, Coverage) %>% matrify() 
Spp_21[] <- lapply(Spp_21, as.numeric)

Spp_22 = dplyr::select(data_22, ID, Species, Coverage) %>% matrify() 
Spp_22[] <- lapply(Spp_22, as.numeric)

# Calculates species richness by year #
SR_20 = specnumber(Spp_20)
SR_20 = as.data.frame(SR_20)

SR_21 = specnumber(Spp_21)
SR_21 = as.data.frame(SR_21)

SR_22 = specnumber(Spp_22)
SR_22 = as.data.frame(SR_22)

# Create Grouped Treatment Table and Summaries to fit Species Table #
Treat_20 <- select(data_20, Treatment, ID) %>% 
  group_by(ID, Treatment) %>% 
  dplyr::summarise()

Treat_21 <- select(data_21, Treatment, ID) %>% 
  group_by(ID, Treatment) %>% 
  dplyr::summarise()

Treat_22 <- select(data_22, Treatment, ID) %>% 
  group_by(ID, Treatment) %>% 
  dplyr::summarise()

## Merge species richness with habitat/plot data for ggplot ##
SR_treat_20 = cbind(Treat_20, SR_20)
colnames(SR_treat_20)[3] <- "SR"
as.numeric(SR_treat_20$SR)
SR_treat_20$SR = as.numeric(SR_treat_20$SR)
SR_treat_20 = as.data.frame(SR_treat_20)

SR_treat_21 = cbind(Treat_21, SR_21)
colnames(SR_treat_21)[3] <- "SR"
as.numeric(SR_treat_21$SR)
SR_treat_21$SR = as.numeric(SR_treat_21$SR)
SR_treat_21 = as.data.frame(SR_treat_21)

SR_treat_22 = cbind(Treat_22, SR_22)
colnames(SR_treat_22)[3] <- "SR"
SR_treat_22$SR = as.numeric(SR_treat_22$SR)
SR_treat_22 = as.data.frame(SR_treat_22)

################################################################################
################ Test for Significance across years ############################
################################################################################

# Check Assumptions #
model  <- lm(SR ~ Treatment, data = SR_treat_20)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
SR_treat_20$Treatment = as.factor(SR_treat_20$Treatment)
SR_treat_20 %>% levene_test(SR_treat_20$SR ~ SR_treat_20$Treatment)

# Test for Significance #
anova_20 = SR_treat_20 %>% kruskal_test(SR ~ Treatment) %>% 
  add_significance()
summary(anova_20)

tukey_20 <- SR_treat_20 %>% 
  dunn_test(SR ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_20

# Check Assumptions #
model  <- lm(SR ~ Treatment, data = SR_treat_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
SR_treat_21$Treatment = as.factor(SR_treat_21$Treatment)
SR_treat_21 %>% levene_test(SR_treat_21$SR ~ SR_treat_21$Treatment)

# Test for Significance #
anova_21 = SR_treat_21 %>% anova_test(SR ~ Treatment) %>% 
  add_significance()
summary(anova_21)

tukey_21 <- SR_treat_21 %>% 
  tukey_hsd(SR ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_21

tmp <- tabular(Treatment ~ SR* (mean+sd+std.error), data=SR_treat_21)
tmp

# Check Assumptions #
model  <- lm(SR ~ Treatment, data = SR_treat_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
SR_treat_22$Treatment = as.factor(SR_treat_22$Treatment)
SR_treat_22 %>% levene_test(SR_treat_22$SR ~ SR_treat_22$Treatment)

# Test for Significance #
anova_22 = SR_treat_22 %>% anova_test(SR ~ Treatment) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- SR_treat_22 %>% 
  tukey_hsd(SR ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

tmp <- tabular(Treatment ~ SR* (mean+sd+std.error), data=SR_treat_22)
tmp

################################################################################
################ Create SR Boxplots for each year ##############################
################################################################################

SR_Box_20 = 
  ggplot(SR_treat_20, aes(x = Treatment, y = SR), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_20,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_20,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_20)) +
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
  labs(x = "", y = "Species Richness", title = "2020")
SR_Box_20

SR_Box_21 = 
  ggplot(SR_treat_21, aes(x = Treatment, y = SR), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_21,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_21,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_21)) +
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
        axis.title.y = ,   
        axis.text.x=element_text(size=15, face = "bold", color = "black"), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Species Richness", title = "2021")
SR_Box_21

SR_Box_22 = 
  ggplot(SR_treat_22, aes(x = Treatment, y = SR), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,
                     hide.ns = T)+
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
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
  labs(x = "", y = "", title = "2022")
SR_Box_22

################## Save Figures Above using ggarrange ##########################
ggarrange(SR_Box_20, SR_Box_21, SR_Box_22, ncol = 2, nrow = 2)
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2020-2022_SR.png", 
       width = 12, height = 7)

