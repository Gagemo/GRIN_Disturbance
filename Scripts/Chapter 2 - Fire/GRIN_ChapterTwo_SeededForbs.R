################################################################################
################################################################################
#########################       GRIN - Fire       ##############################
#########################       Native Forbs      ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2023        ##############################
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

####################### Read in 2021 - 2023 Data  ##############################
GRIN = read.csv("Data/GRIN - 2021-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Select only Seeding Treatment # 
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

data = filter(GRIN, Species == "Pityopsis trayci" | 
              Species == "Liatris gracilis")

# Creates data sets by year #
data_22 = filter(data, Year == 2)
data_23 = filter(data, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Data ######################################

# Sums CV data #
GRIN_22M <- dplyr::select(data_22, ID, Group, Coverage) %>% matrify()

# Create grouped treatment/environment table and summaries to fit species table#
Treat_22 = group_by(data_22, ID, Fire) %>% dplyr::summarize()
Forb_22 = select(GRIN_22M, Forb) %>%
  cbind(Treat_22)

# Check Assumptions #
model  <- lm(Forb ~ Fire, data = Forb_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb_22$Fire = as.factor(Forb_22$Fire)
Forb_22 %>% levene_test(Forb ~ Fire)

# Test for Significance #
anova_22 = Forb_22 %>% kruskal_test(Forb ~ Fire) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Forb_22 %>% 
  dunn_test(Forb ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Data ######################################
# Sums CV data #
GRIN_23M <- dplyr::select(data_23, ID, Group, Coverage) %>% matrify()

# Create grouped treatment/environment table and summaries to fit species table#
Treat_23 = group_by(data_23, ID, Fire) %>% dplyr::summarize()

Forb_23 = select(GRIN_23M, Forb) %>%
  cbind(Treat_23)

# Check Assumptions #
model  <- lm(Forb ~ Fire, data = Forb_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Forb_23$Fire = as.factor(Forb_23$Fire)
Forb_23 %>% levene_test(Forb ~ Fire)

# Test for Significance #
anova_23 = Forb_23 %>% kruskal_test(Forb ~ Fire) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Forb_23 %>% 
  dunn_test(Forb ~ Fire) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

################################################################################
################ Create Box Plot for Seeded forbs Across Years #################
################################################################################

## Lovegrass Coverage 2022 Box plot ##
Box22 = 
  ggplot(Forb_22, aes(x = Fire, y = Forb), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,
                     hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_22,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_22)) +
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
  labs(x = "", y = "Seeded Forb % Coverage", title = "2022")
    Box22

ggsave("Figures/Chapter 2 - Fire/seedbox22.png", 
       width = 12, height = 8)

## Lovegrass Coverage 2023 Boxplot ##
Box23 = 
  ggplot(Forb_23, aes(x = Fire, y = Forb), colour = Fire) +
  geom_boxplot(aes(fill=Fire), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Fire), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_23,
                     hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
        axis.title.y = element_blank(),   
        axis.text.x=element_text(size=15, face = "bold", color = "black"),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.ticks = element_blank(), 
        strip.text.x = 
          element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "", title = "2023")
Box23

ggsave("Figures/Chapter 2 - Fire/seedbox23.png", 
       width = 12, height = 8)

################## Save Figures Above using ggarrange ##########################
ggarrange(Box22, Box23, ncol = 2, nrow = 1)
ggsave("Figures/Chapter 2 - Fire/22-23_seedBox.png", 
       width = 12, height = 10)
