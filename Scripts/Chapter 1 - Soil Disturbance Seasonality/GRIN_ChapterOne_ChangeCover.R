################################################################################
################################################################################
#########################    Data - Disturbance   ##############################
#########################    Change in Cover      ##############################
#########################  University of Florida  ##############################
#########################     Gage LaPierre       ##############################
#########################      2021 - 2022        ##############################
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
                         by = c('ID_', "Treatment", 'Species'))
Abundance_w = arrange(Abundance_w, Treatment)

# Turns NA values into zeros #
Abundance_w$Coverage.x <- ifelse(is.na(Abundance_w$Coverage.x), 0, 
                                 Abundance_w$Coverage.x)
Abundance_w$Coverage.y <- ifelse(is.na(Abundance_w$Coverage.y), 0, 
                                 Abundance_w$Coverage.y)

# Change abundance to reflect percentage change from (Year 1) to (Year 2)  #
Change_Abundance <- Abundance_w %>% 
  dplyr::select(ID_, Treatment, Species, 
                Coverage.x, Coverage.y) %>%
  group_by(ID_, Treatment, Species) %>% 
  mutate(Change_abundance = Coverage.y - Coverage.x)

##################################  COVER CAHNGIh ##############################

################################################################################
########################### Hairy Indig ########################################
################################################################################
Ih = 
  Change_Abundance[which(Change_Abundance$Species == "Indigofera hirsuta"),]
Ih<-as.data.frame(Ih)
Ih$Treatment<-factor(Ih$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Ih)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Ih = Ih %>% anova_test(Change_abundance ~ Treatment) %>% 
  add_significance()
anova_Ih

lm(formula = Change_abundance ~ Treatment, Ih)
tukey_Ih <- Ih %>% 
  tukey_hsd(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Ih

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Ih)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ih_Change.csv")

Ih_Change_Box = 
  ggplot(Ih, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Ih,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Ih, detailed = TRUE),
       caption = get_pwc_label(tukey_Ih)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Indigo hirsuta")
Ih_Change_Box

################################################################################
########################### Oxalis corniculata #################################
################################################################################
Oc = 
  Change_Abundance[which(Change_Abundance$Species == "Oxalis corniculata"),]
Oc<-as.data.frame(Oc)
Oc$Treatment<-factor(Oc$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Oc)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Oc = Oc %>% anova_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Oc)

tukey_Oc <- Oc %>% 
  tukey_hsd(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Oc

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Oc )
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Oc_Change.csv")

Oc_change_Box = 
  ggplot(Oc, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Oc,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Oc, detailed = TRUE),
       caption = get_pwc_label(tukey_Oc)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Oxalis corniculata")
Oc_change_Box

 ################################################################################
########################### Paspalum setaceum ##########################################
################################################################################
Ps = 
  Change_Abundance[which(Change_Abundance$Species == "Paspalum setaceum"),]
Ps<-as.data.frame(Ps)
Ps$Treatment<-factor(Ps$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Ps)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Ps = Ps %>% anova_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Ps)

tukey_Ps <- Ps %>% 
  tukey_hsd(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Ps

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Ps)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ps_Change.csv")

Ps_change_Box = 
  ggplot(Ps, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Ps,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Ps, detailed = TRUE),
       caption = get_pwc_label(tukey_Ps)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Paspalum setaceum")
Ps_change_Box

################################################################################
############################### Rubus ##########################################
################################################################################
Rubus = 
  Change_Abundance[which(Change_Abundance$Species == "Rubus spp."),]
Rubus<-as.data.frame(Rubus)
Rubus$Treatment<-factor(Rubus$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Rubus)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Rubus = Rubus %>% kruskal_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Rubus)

tukey_Rubus <- Rubus %>% 
  dunn_test(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Rubus

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Rubus)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Rubus_Change.csv")

Rubus_change_Box = 
  ggplot(Rubus, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Rubus,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Rubus, detailed = TRUE),
       caption = get_pwc_label(tukey_Rubus)) +
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
  labs(x = "", y = "% Change in Coverage", title = "Rubus spp.")
Rubus_change_Box

################################################################################
############################### Paspalum notatum ########################################
################################################################################
Pn = 
  Change_Abundance[which(Change_Abundance$Species == "Paspalum notatum"),]
Pn<-as.data.frame(Pn)
Pn$Treatment<-factor(Pn$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Pn)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Pn = Pn %>% kruskal_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Pn)

tukey_Pn <- Pn %>% 
  dunn_test(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Pn
summary(tukey_Pn)

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Pn)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_Change.csv")

Pn_change_Box = 
  ggplot(Pn, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Pn,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Pn, detailed = TRUE),
       caption = get_pwc_label(tukey_Pn)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Paspalum notatum")
Pn_change_Box
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/change_Pn.png", 
       width = 10, height = 7)

################################################################################
############################### Richardia spp.   ########################################
################################################################################
Rs = 
  Change_Abundance[which(Change_Abundance$Species == "Richardia spp."),]
Rs<-as.data.frame(Rs)
Rs$Treatment<-factor(Rs$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Rs)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Rs = Rs %>% kruskal_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Rs)

tukey_Rs <- Rs %>% 
  dunn_test(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Rs
summary(tukey_Rs)

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Rs)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Rs_Change.csv")

Rs_change_Box = 
  ggplot(Rs, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Rs,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Rs, detailed = TRUE),
       caption = get_pwc_label(tukey_Rs)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Richardia spp.")
Rs_change_Box

################################################################################
############################### Paspalum notatum ########################################
################################################################################
Pn = 
  Change_Abundance[which(Change_Abundance$Species == "Paspalum notatum"),]
Pn<-as.data.frame(Pn)
Pn$Treatment<-factor(Pn$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Pn)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Pn = Pn %>% kruskal_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Pn)

tukey_Pn <- Pn %>% 
  dunn_test(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Pn
summary(tukey_Pn)

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Pn)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Pn_Change.csv")

Pn_change_Box = 
  ggplot(Pn, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Pn,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Pn, detailed = TRUE),
       caption = get_pwc_label(tukey_Pn)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Paspalum notatum")
Pn_change_Box
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/change_Pn.png", 
       width = 10, height = 7)

################################################################################
############################### Dichondra spp.   ###############################
################################################################################
Ds = 
  Change_Abundance[which(Change_Abundance$Species == "Dichondra spp."),]
Ds<-as.data.frame(Ds)
Ds$Treatment<-factor(Ds$Treatment)

# Check Assumptions #
model  <- lm(Change_abundance ~ Treatment, data = Ds)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
plot(model, 1)

# Test for Significance #
anova_Ds = Ds %>% kruskal_test(Change_abundance ~ Treatment) %>% 
  add_significance()
summary(anova_Ds)

tukey_Ds <- Ds %>% 
  dunn_test(Change_abundance ~ Treatment) %>% 
  add_significance() %>% 
  add_xy_position()
tukey_Ds
summary(tukey_Ds)

tmp <- tabular(Treatment ~ Change_abundance* (mean+sd+std.error), data=Ds)
tmp

write.csv.tabular(tmp, "Figures/Chapter 1 - Soil Disturbance Seasonality/Ds_Change.csv")

Ds_change_Box = 
  ggplot(Ds, aes(x = Treatment, y = Change_abundance), colour = Treatment) +
  geom_boxplot(aes(fill=Treatment), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(fill=Treatment), size = 3, 
             position = position_jitterdodge(), alpha = 0.7) +
  stat_pvalue_manual(tukey_Ds,size = 8, bracket.size = 1, hide.ns = T)+
  labs(subtitle = get_test_label(anova_Ds, detailed = TRUE),
       caption = get_pwc_label(tukey_Ds)) +
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
  labs(x = "", y = "Change in % Coverage", title = "Dichondra spp.")
Ds_change_Box


################## Save Figures Above using ggarrange ##########################
Change = 
  ggarrange(Ih_Change_Box, Oc_change_Box, Rs_change_Box, 
            Ps_change_Box, Pn_change_Box, Ds_change_Box, ncol = 2, nrow = 3)
Change
annotate_figure(Change, top = text_grob("", color = "black", 
                                        face = "bold", size = 25))
ggsave("Figures/Chapter 1 - Soil Disturbance Seasonality/Change.png", 
       width = 12, height = 14)

