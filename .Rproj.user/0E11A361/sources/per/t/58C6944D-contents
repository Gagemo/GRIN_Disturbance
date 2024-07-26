################################################################################
################################################################################
#########################   GRIN - Disturbance    ##############################
#########################       Bare Ground        #############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2022       ###############################
################################################################################
################################################################################

######################### Clears Environment & History  ########################
rm(list=ls(all=TRUE))
cat("\014") 

#########################     Installs Packages   ##############################
list.of.packages <- c("tidyverse", "vegan", "agricolae", "tables", "plotrix")
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

##########################     Loads Packages     ##############################
library(tidyverse)
library(vegan)
library(agricolae)
library(tables)
library(plotrix)

##########################     Read in 2022 Data  ##############################
GRIN = read.csv("Data/GRIN - 2020-2023.csv")
GRIN$Coverage = as.numeric(GRIN$Coverage)
GRIN$Plot = as.character(GRIN$Plot)

str(GRIN)
summary(GRIN)

# Remove Seeding Treatment # 
GRIN = filter(GRIN, Treatment != 'S')

# Remove Year Three #
GRIN = filter(GRIN, Year != 3)

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
  grepl(9, Coverage) ~ 85
))

# Filters data for just bare ground #
bare = filter(GRIN, Group == "Bare")

# Creates data sets by year #
bare_21 = filter(bare, Year == 1)
bare_22 = filter(bare, Year == 2)

################################################################################
################ Test for Significance across years ############################
############################### 2021 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = bare_21)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
bare_21$Treatment = as.factor(bare_21$Treatment)
bare_21 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_bare21 = bare_21 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_bare21)

tukey_bare21 <- bare_21 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_bare21

############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ Treatment, data = bare_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
bare_22$Treatment = as.factor(bare_22$Treatment)
bare_22 %>% levene_test(Coverage ~ Treatment)

# Test for Significance #
anova_bare22 = bare_22 %>% kruskal_test(Coverage ~ Treatment) %>% 
  add_significance()
summary(anova_bare22)

tukey_bare22 <- bare_22 %>% 
  dunn_test(Coverage ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_bare22

## Bare Ground Coverage 2021 ##
BareBox21 = 
  ggplot(bare_21, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 4, alpha = 0.5) +
  stat_pvalue_manual(tukey_bare21, hide.ns = T, size = 8, bracket.size = 1)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_bare21, detailed = TRUE),
       caption = get_pwc_label(tukey_bare21)) +
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
  labs(x = "Treatment", y = "Bare Ground % Coverage", title = "2021")
BareBox21

ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/Bareground_box21.png", 
       width = 10, height = 7)

## Bare Ground Coverage 2022 ##
BareBox22 = 
  ggplot(bare_22, aes(x = Treatment, y = Coverage), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), 
             position = position_jitterdodge(), size = 4, alpha = 0.5) +
  stat_pvalue_manual(tukey_bare22, hide.ns = T, size = 8, bracket.size = 1) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_bare22, detailed = TRUE),
       caption = get_pwc_label(tukey_bare22)) +
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
  labs(x = "Treatment", y = "", title = "2022")
BareBox22

ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/Bareground_box22.png", 
       width = 10, height = 7)

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=bare_21)
tmp

tmp <- tabular(Treatment ~ Coverage* (mean+sd+std.error), data=bare_22)
tmp

################## Save Figures Above using ggarrange ##########################
ggarrange(BareBox21, BareBox22, ncol = 2, nrow = 1)
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/2021-2022_Bare.png", 
       width = 12, height = 7)

