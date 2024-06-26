################################################################################
################################################################################
#########################       GRIN - Fire        #############################
#########################       Bare Ground        #############################
#########################  University of Florida   #############################
#########################     Gage LaPierre        #############################
#########################      2021 - 2022         #############################
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

##########################     Read in 2022 Data  ##############################
GRIN = read.csv("Data/GRIN - 2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Select just Seeding Treatment # 
GRIN = filter(GRIN, Treatment == 'S')

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

# Filters data for just bare ground #
bare = filter(GRIN, Group == "Bare")

# Creates data sets by year #
bare_22 = filter(bare, Year == 2)
bare_23 = filter(bare, Year == 3)

################################################################################
################ Test for Significance across years ############################
############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Fire, data = bare_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
bare_22$Fire = as.factor(bare_22$Fire)
bare_22 %>% levene_test(Coverage ~ Fire)

# Test for Significance #
anova_bare22 = bare_22 %>% kruskal_test(Coverage ~ Fire) %>% 
  add_significance()
summary(anova_bare22)

tukey_bare22 <- bare_22 %>% 
  dunn_test(Coverage ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_bare22

############################### 2023 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Fire, data = bare_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
bare_23$Fire = as.factor(bare_23$Fire)
bare_23 %>% levene_test(Coverage ~ Fire)

# Test for Significance #
anova_bare23 = bare_23 %>% kruskal_test(Coverage ~ Fire) %>% 
  add_significance()
summary(anova_bare23)

tukey_bare23 <- bare_23 %>% 
  dunn_test(Coverage ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_bare23

## Bare Ground Coverage 2022 ##
BareBox22 = 
  ggplot(bare_22, aes(x = Fire, y = Coverage), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_bare22,
                     hide.ns = T)+
  ylim(0, 70) +
  labs(subtitle = get_test_label(anova_bare22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_bare22)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Till', 'Late-Spring', 'Winter')) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", colour = "black"),
        text=element_text(size=16),
        axis.title.x = element_text(size=15, face="bold", colour = "black"),    
        axis.title.y = element_text(size=15, face="italic", colour = "black"),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_text(size=15, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Bare Ground % Coverage", title = "2022")
BareBox22

ggsave("Figures/Chapter 2 - Fire/Bareground_box22.png", 
       width = 12, height = 8)

## Bare Ground Coverage 2023 ##
BareBox23 = 
  ggplot(bare_23, aes(x = Fire, y = Coverage), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), position = position_jitterdodge(), 
             size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_bare23, hide.ns = T, size = 8)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_bare23,detailed = TRUE),
       caption = get_pwc_label(tukey_bare23)) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
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
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "Bare Ground % Coverage", title = "2023")
BareBox23

ggsave("Figures/Chapter 2 - Fire/Bareground_box23.png", 
       width = 12, height = 8)

################## Save Figures Above using ggarrange ##########################
ggarrange(BareBox22, BareBox23, ncol = 2, nrow = 1)
ggsave("Figures/Chapter 2 - Fire/22-23_Bare.png", 
       width = 12, height = 10)

