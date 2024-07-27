################################################################################
################################################################################
#########################      GRIN - Fire        ##############################
#########################   Species Coverage      ##############################
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

# Select only Seeding Treatment # 
Data = filter(Data, Treatment == 'S')

################################################################################
################################################################################
################################################################################
############################### Pityopsis trayci ###############################
################################################################################
################################################################################
################################################################################

Pt = filter(Data, Species == "Pityopsis trayci")
summary(Pt)

# Creates data sets by year #
Pt_22 = filter(Pt, Year == 2)
Pt_23 = filter(Pt, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Pt_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Pt_22$T.F= as.factor(Pt_22$T.F)
Pt_22 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_22 = Pt_22 %>% anova_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Pt_22 %>% 
  tukey_hsd(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Pt ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Pt_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Pt_23$T.F= as.factor(Pt_23$T.F)
Pt_23 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_23 = Pt_23 %>% anova_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Pt_23 %>% 
  tukey_hsd(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

## Lovegrass Coverage 2022 Box plot ##
PtBox22 = 
  ggplot(Pt_22, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
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
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "P. trayci % Coverage", title = "2022")
PtBox22

## Pityopsis trayci  Coverage 2023 Boxplot ##
PtBox23 = 
  ggplot(Pt_23, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "P. trayci % Coverage", title = "2023")
PtBox23

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Pt_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Pt_22.csv")

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Pt_23)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Pt_23.csv")

################################################################################
################################################################################
############################### Liatris ########################################
################################################################################
################################################################################

Lg = filter(Data, Species == "Liatris gracilis")
summary(Lg)

# Creates data sets by year #
Lg_22 = filter(Lg, Year == 2)
Lg_23 = filter(Lg, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Data ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Lg_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Lg_22$T.F= as.factor(Lg_22$T.F)
Lg_22 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_22 = Lg_22 %>% anova_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Lg_22 %>% 
  tukey_hsd(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Lg ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Lg_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Lg_23$T.F= as.factor(Lg_23$T.F)
Lg_23 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_23 = Lg_23 %>% anova_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Lg_23 %>% 
  tukey_hsd(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

##  Coverage 2022 Box plot ##
LgBox22 = 
  ggplot(Lg_22, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "L. gracilis % Coverage", title = "2022")
LgBox22

## Coverage 2023 Boxplot ##
LgBox23 = 
  ggplot(Lg_23, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "L. gracilis % Coverage", title = "2023")
LgBox23

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Lg_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Lg_22.csv")

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Lg_23)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Lg_23.csv")

################################################################################
################################################################################
###################### Sorghastrum #############################################
################################################################################
################################################################################

Ss = filter(Data, Species == "Sorghastrum secundum")
summary(Ss)

# Creates data sets by year #
Ss_22 = filter(Ss, Year == 2)
Ss_23 = filter(Ss, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Ss ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Ss_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ss_22$T.F= as.factor(Ss_22$T.F)
Ss_22 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_22 = Ss_22 %>% anova_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Ss_22 %>% 
  tukey_hsd(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Ss ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Ss_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ss_23$T.F= as.factor(Ss_23$T.F)
Ss_23 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_23 = Ss_23 %>% anova_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Ss_23 %>% 
  tukey_hsd(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

## indiangrass Coverage 2022 Box plot ##
SsBox22 = 
  ggplot(Ss_22, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
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
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "S. secundum % Coverage", title = "2022")
SsBox22

## Coverage 2023 Boxplot ##
SsBox23 = 
  ggplot(Ss_23, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "S. secundum % Coverage", title = "2023")
SsBox23

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Ss_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Ss_22.csv")

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Ss_23)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Ss_23.csv")

################################################################################
################################################################################
################################ Splitbeard ####################################
################################################################################
################################################################################

At = filter(Data, Species == "Andropogon ternarius")
summary(At)

# Creates data sets by year #
At_22 = filter(At, Year == 2)
At_23 = filter(At, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 At ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = At_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
At_22$T.F= as.factor(At_22$T.F)
At_22 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_22 = At_22 %>% kruskal_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- At_22 %>% 
  dunn_test(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 At ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = At_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
At_23$T.F= as.factor(At_23$T.F)
At_23 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_23 = At_23 %>% kruskal_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- At_23 %>% 
  dunn_test(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

## Coverage 2022 Box plot ##
AtBox22 = 
  ggplot(At_22, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "A. ternarius % Coverage", title = "2022")
AtBox22

ggsave("Figures/Chapter 2 - Fire/Split_box22.png", 
       width = 12, height = 8)

## Coverage 2023 Boxplot ##
AtBox23 = 
  ggplot(At_23, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "A. ternarius % Coverage", title = "2023")
AtBox23

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=At_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/At_22.csv")

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=At_23)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/At_23.csv")

################################################################################
################################################################################
################################ Wiregrass  ####################################
################################################################################
################################################################################

Ab = filter(Data, Species == "Aristida beyrichiana")
summary(Ab)

# Creates data sets by year #
Ab_22 = filter(Ab, Year == 2)
Ab_23 = filter(Ab, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Ab ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Ab_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ab_22$T.F= as.factor(Ab_22$T.F)
Ab_22 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_22 = Ab_22 %>% kruskal_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Ab_22 %>% 
  dunn_test(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Ab ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Ab_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Ab_23$T.F= as.factor(Ab_23$T.F)
Ab_23 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_23 = Ab_23 %>% kruskal_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Ab_23 %>% 
  dunn_test(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

## Coverage 2022 Box plot ##
AbBox22 = 
  ggplot(Ab_22, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "A. beyrichiana % Coverage", title = "2022")
AbBox22

## Coverage 2023 Boxplot ##
AbBox23 = 
  ggplot(Ab_23, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "A. beyrichiana % Coverage", title = "2023")
AbBox23

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Ab_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Ab_22.csv")

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Ab_23)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Ab_23.csv")

################################################################################
################################################################################
################################ Lovegrass  ####################################
################################################################################
################################################################################

Es = filter(Data, Species == "Eragrostis spectabilis")
summary(Es)

# Creates data sets by year #
Es_22 = filter(Es, Year == 2)
Es_23 = filter(Es, Year == 3)

################################################################################
################ Test for Significance across years ############################
################################################################################

############################### 2022 Es ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Es_22)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Es_22$T.F= as.factor(Es_22$T.F)
Es_22 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_22 = Es_22 %>% kruskal_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_22)

tukey_22 <- Es_22 %>% 
  dunn_test(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_22

############################### 2023 Es ######################################
# Check Assumptions #
model  <- lm(Coverage ~ T.F, data = Es_23)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)
# Compute Levene's Test
Es_23$T.F= as.factor(Es_23$T.F)
Es_23 %>% levene_test(Coverage ~ T.F)

# Test for Significance #
anova_23 = Es_23 %>% kruskal_test(Coverage ~ T.F) %>% 
  add_significance()
summary(anova_23)

tukey_23 <- Es_23 %>% 
  dunn_test(Coverage ~ T.F) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_23

## Coverage 2022 Box plot ##
EsBox22 = 
  ggplot(Es_22, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_22,size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
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
        axis.title.y = element_text(size=20, face="italic", colour = "black"),   
        axis.text.x=element_text(size=20, face = "bold", color = "black"),
        axis.text.y=element_text(size=20, face = "bold", color = "black"),
        strip.text.x = 
          element_text(size = 20, colour = "black", face = "bold"),
        legend.position = "none") +
  guides(fill = guide_legend(label.position = "bottom")) +
  labs(x = "", y = "E. spectabilis % Coverage", title = "2022")
EsBox22

## Coverage 2023 Boxplot ##
EsBox23 = 
  ggplot(Es_23, aes(x = T.F, y = Coverage), colour = T.F) +
  geom_boxplot(aes(fill=T.F), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=T.F), 
             position = position_jitterdodge(), size = 2, alpha = 0.5) +
  stat_pvalue_manual(tukey_23, size = 8, bracket.size = 1, hide.ns = T)+
  scale_fill_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                    values=c("#333333", "#FF9900", "#3366FF")) +
  scale_color_manual(labels=c('No Burn', 'Late-Spring', 'Winter'),
                     values=c("#333333", "#FF9900", "#3366FF")) +
  scale_x_discrete(labels=c('No-Burn', 'Late-Spring', 'Winter')) +
  ylim(0, 100) +
  labs(subtitle = get_test_label(anova_23,
                                 detailed = TRUE),
       caption = get_pwc_label(tukey_23)) +
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
  labs(x = "", y = "E. spectabilis % Coverage", title = "2023")
EsBox23

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Es_22)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Es_22.csv")

tmp <- tabular(T.F ~ Coverage * (mean+sd+std.error), data=Es_23)
tmp

write.csv.tabular(tmp, "Figures/Chapter 2 - Fire/Es_23.csv")

################## Save All Figures/Chapter 2 - Fire Above using ggarrange ##########################
all = 
  ggarrange(PtBox22, PtBox23, LgBox22, LgBox23, EsBox22, EsBox23,
            AbBox22, AbBox23, SsBox22, SsBox23, AtBox22, AtBox23,
             ncol = 2, nrow = 6)
all
ggsave("Figures/Chapter 2 - Fire/22-23_Obli_Forb_Grass.png", 
       width = 14, height = 22)

