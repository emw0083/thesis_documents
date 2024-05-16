### Identifying mean egg diameter for the three dams used ####


#Packages
library(readxl)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
cbPalette <- c("#643d6e", "#a982b4", "#d4c0d9")
library(tidyverse)
library(rstatix)
library(lme4) 
library(ggpubr)
library(afex)
library("gridExtra")
library(emmeans)


####------------------General Data Management----------------####

#Read in file
raw <- read_excel("Egg_measurements.xlsx")

#makes new df with columns needed 
data <- raw[ , c("Tube ID", "Sample ID", "Photo ID", "Diameter (um)" )]
#View(data)

#Renames columns with no spaces
names(data)[names(data) == 'Tube ID'] <- 'tube'
names(data)[names(data) == 'Sample ID'] <- 'ID'
names(data)[names(data) == 'Photo ID'] <- 'sampleID'
names(data)[names(data) == 'Diameter (um)'] <- 'diameter'

#Omits data with NA
df = na.omit(data)
#View(df)

#Takes the mean of the three measurements for each egg sampled 
means = aggregate(df$diameter, list(df$sampleID, df$ID, df$tube), mean) 

#Renames columns to reflect what they are 
colnames(means)[colnames(means) == "x"] = "avg_diameter"
colnames(means)[colnames(means) == "Group.1"] = "Individual"
colnames(means)[colnames(means) == "Group.2"] = "ID"
colnames(means)[colnames(means) == "Group.3"] = "tube"

#Adds column with this information
means$Dam <- substr(means$tube,6,6)
means$Replicate <- substr(means$ID,6,8)

#file for comparing pluteus length and pigment cell count to egg diameter
#write_csv(means, "mean_egg_diameter.csv")


#Renames 1, 2, 3 to Dam 1, Dam 2, and Dam 3 
means["Dam"][means["Dam"]=="1"]<- "Dam 1"
means["Dam"][means["Dam"]=="2"]<- "Dam 2"
means["Dam"][means["Dam"]=="3"]<- "Dam 3"

####-------------STATS---------------------####

lmer_lengthdam=lmer(avg_diameter ~ Dam + (1|ID), data= means)
summary(lmer_lengthdam)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: avg_diameter ~ Dam + (1 | ID)
#    Data: means
# 
# REML criterion at convergence: 1554.4
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.2905 -0.6213  0.0265  0.6041  3.4544 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  ID       (Intercept)  1.601   1.265   
#  Residual             19.954   4.467   
# Number of obs: 266, groups:  ID, 9
# 
# Fixed effects:
#             Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)  94.1726     0.8692  5.8133 108.343 7.83e-11 ***
# DamDam 2     -5.3101     1.2337  5.8961  -4.304  0.00528 ** 
# DamDam 3     -6.0349     1.2292  5.8133  -4.909  0.00293 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#          (Intr) DamDm2
# DamDam 2 -0.705       
# DamDam 3 -0.707  0.498
# 
# emm = emmeans(lmer_lengthfem, specs = pairwise ~ Dam)
# summary(emm)
# #
# 
# 
# # $emmeans
# # Dam   emmean    SE   df lower.CL upper.CL
# # Dam 1   94.2 0.869 5.94     92.0     96.3
# # Dam 2   88.9 0.876 6.11     86.7     91.0
# # Dam 3   88.1 0.869 5.94     86.0     90.3
# # 
# # Degrees-of-freedom method: kenward-roger 
# # Confidence level used: 0.95 
# # 
# # $contrasts
# # contrast            estimate   SE   df t.ratio p.value
# # Dam 1 - Dam 2    5.310 1.23 6.03   4.304  0.0119
# # Dam 1 - Dam 3    6.035 1.23 5.94   4.909  0.0066
# # Dam 2 - Dam 3    0.725 1.23 6.03   0.587  0.8317
# # 
# # Degrees-of-freedom method: kenward-roger 
# # P value adjustment: tukey method for comparing a family of 3 estimates 


####------------------Figures----------------####
p1 = ggplot(means, aes(x=Dam, y=avg_diameter, fill=Dam)) +
  geom_boxplot(width=0.6) +
  geom_point(position=position_jitterdodge(0.25), shape=21, size=3) +
  scale_fill_manual(values=cbPalette) +
  ylab("Mean Egg Diameter (micron)") + 
  theme(legend.position="none",
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title=element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.75)) + scale_y_continuous(limits=c(70,115)) +
  geom_signif(comparisons=list(c("Dam 1", "Dam 2")), annotations= "**", y_position = 105.75, textsize=9) +
  geom_signif(comparisons=list(c("Dam 1", "Dam 3")), annotations= "**", y_position = 109.75, textsize=9) 

p1 

p1 + theme(plot.title = element_text(face = "bold"))


ggsave("egg_plot_updated_pretty.svg", p1)




##-------------Stats Table--------------

library(sjstats)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(parameters)
library(sjlabelled)


tab_model(lmer_lengthdam)
tab_model(emm)


library(report)
report(lmer_lengthdam)

# We fitted a linear mixed model (estimated using REML and nloptwrap optimizer) to predict
# avg_diameter with Dam (formula: avg_diameter ~ Dam). The model included ID as random
# effect (formula: ~1 | ID). The model's total explanatory power is substantial
# (conditional R2 = 0.31) and the part related to the fixed effects alone (marginal R2) is
# of 0.25. The model's intercept, corresponding to Dam = Dam 1, is at 94.17 (95% CI [92.46, 95.88], t(261) = 108.34, p < .001). Within this model:
# - The effect of Dam [Dam 2] is statistically significant and negative (beta = -5.31, 95%
#                                                                          CI [-7.74, -2.88], t(261) = -4.30, p < .001; Std. beta = -1.00, 95% CI [-1.46, -0.54])
# - The effect of Dam [Dam 3] is statistically significant and negative (beta = -6.03, 95%
#                                                                        CI [-8.46, -3.61], t(261) = -4.91, p < .001; Std. beta = -1.14, 95% CI [-1.59, -0.68])
# 
# Standardized parameters were obtained by fitting the model on a standardized version of
# the dataset. 95% Confidence Intervals (CIs) and p-values were computed using a Wald
# t-distribution approximation.





