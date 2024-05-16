#######----- Consolidated pigment cell data for manuscript----- #############


#Load Packages
library(readxl)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
cbPalette <- c("#0072B2","darkgrey","#84B9EF","#B73E3E", "#DD5353", "#D61C4E", "#FFACAC", "#952E4B", "#0072B2","#CC79A7", "#56B4E9", "#E69F00","#999999", "#009E73", "#F0E442", "#D55E00", "#CC79A7")
cbPalette1 <- c("#84B9EF","#B73E3E")
library(tidyverse)
library(rstatix)
library(lme4) #Needed for stats
library(wesanderson)
library(ggpubr)
library(afex) #Needed for stats! (this gives the p value)
library("gridExtra")
library(cowplot)
library(dplyr)
library(plyr)
library(multcomp) ####lmer post hoc help
library(emmeans)

#for stats table
library(sjstats)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(parameters)
library(sjlabelled)

####-------------General File manipulation-------------####

##Read in datafile
raw <- read_excel("Pigment_Cell_and_Length_Data.xlsx")
#View(raw)

##Subset the data so only working with values needed for pigment cell analysis
data <- raw[ , c("Tube_ID","Sample_ID", "Photo_ID", "Counter_Type", "Count")]
#view(data)

##Subset further so it's just Type 1s, which is the label for the pigment cell point
types <- subset(data, Counter_Type == "Type1", select = c("Tube_ID", "Sample_ID", "Photo_ID", "Counter_Type", "Count"))
#view(types)

##Renames columns so they have shorter names to work with
names(types)[names(types) == 'Tube_ID'] <- 'tube'
names(types)[names(types) == 'Photo_ID'] <- 'stack'
names(types)[names(types) == 'Counter_Type'] <- 'type'
names(types)[names(types) == 'Count'] <- 'count'


##Omits rows with NA so they are removed from analysis
df = na.omit(types)
#view(df)

##Need to write csv of pigment cell subset for future correlation figure
#write.csv(df, "Pigment_subset_CORRECT.csv")


##Extracts values from columns and makes new columns for easier analysis for tested conditions

df$Treatment <- substr(df$tube,2,2)
df$Culture <- substr(df$tube,3,5)
df$Dam <- substr(df$Culture,1,1)
df$Sire <- substr(df$Culture,2,2)
df$Combo <-substr(df$Culture,1,2)

##Rename E and A to Elevated (18°C) and Ambient (14°C)
df["Treatment"][df["Treatment"]=="E"]<- "Elevated (18°C)"
df["Treatment"][df["Treatment"]=="A"]<- "Ambient (14°C)"
#view(df)

#write.csv(df, "Pigment_subset_CORRECT.csv")

#data = read_csv("Pigment_subset_CORRECT.csv")
#View(data)


#####----------------------------STATS-------------------------#####
library(lme4) 

#For Combo --- do not use this, need to use interaction (based on anova test)
lmer_pig=lmer(count ~ Treatment+Combo + (1|Culture:Treatment), data= df)
summary(lmer_pig)

##USE THIS ONE
lmer_pig1=lmer(count ~ Treatment + Combo + Treatment*Combo + (1|Culture:Treatment), data= df)
summary(lmer_pig1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: count ~ Treatment + Combo + Treatment * Combo + (1 | Culture:Treatment)
#    Data: df
# 
# REML criterion at convergence: 2907.2
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.5514 -0.6396 -0.1054  0.5714  4.2492 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 11.23    3.351   
#  Residual                      84.38    9.186   
# Number of obs: 402, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                                  Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                        52.037      2.664   6.810  19.532 3.11e-07 ***
# TreatmentElevated (18°C)           23.113      3.850   7.425   6.004 0.000431 ***
# Combo12                             2.093      3.761   6.760   0.556 0.595854    
# Combo23                            -5.603      3.757   6.737  -1.491 0.181161    
# Combo33                            -6.503      3.757   6.737  -1.731 0.128776    
# TreatmentElevated (18°C):Combo12   -7.489      5.463   7.522  -1.371 0.209954    
# TreatmentElevated (18°C):Combo23  -15.597      5.372   7.040  -2.903 0.022743 *  
# TreatmentElevated (18°C):Combo33  -18.313      5.990   6.964  -3.057 0.018515 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Comb12 Comb23 Comb33 TE(18°C):C1 TE(18°C):C2
# TrtmE(18°C) -0.692                                                       
# Combo12     -0.708  0.490                                                
# Combo23     -0.709  0.491     0.502                                      
# Combo33     -0.709  0.491     0.502  0.503                               
# TE(18°C):C1  0.488 -0.705    -0.688 -0.346 -0.346                        
# TE(18°C):C2  0.496 -0.717    -0.351 -0.699 -0.352  0.505                 
# TE(18°C):C3  0.445 -0.643    -0.315 -0.315 -0.627  0.453       0.461    


anova(lmer_pig1, lmer_pig)

# Data: df
# Models:
#   lmer_pig: count ~ Treatment + Combo + (1 | Culture:Treatment)
# lmer_pig1: count ~ Treatment + Combo + Treatment * Combo + (1 | Culture:Treatment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# lmer_pig     7 2964.5 2992.5 -1475.2   2950.5                        
# lmer_pig1   10 2955.0 2995.0 -1467.5   2935.0 15.458  3   0.001465 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


em_combo = emmeans(lmer_pig1, specs = pairwise ~ Combo|Treatment)
summary(em_combo)

# $emmeans
# Treatment = Ambient (14°C):
#   Combo emmean   SE   df lower.CL upper.CL
# 11      52.0 2.66 6.71     45.7     58.4
# 12      54.1 2.65 6.62     47.8     60.5
# 23      46.4 2.65 6.57     40.1     52.8
# 33      45.5 2.65 6.57     39.2     51.9
# 
# Treatment = Elevated (18°C):
#   Combo emmean   SE   df lower.CL upper.CL
# 11      75.2 2.78 7.95     68.7     81.6
# 12      69.8 2.83 8.47     63.3     76.2
# 23      54.0 2.65 6.57     47.6     60.3
# 33      50.3 3.75 6.57     41.4     59.3
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast          estimate   SE   df t.ratio p.value
# Combo11 - Combo12    -2.09 3.76 6.66  -0.556  0.9417
# Combo11 - Combo23     5.60 3.76 6.64   1.491  0.4918
# Combo11 - Combo33     6.50 3.76 6.64   1.731  0.3798
# Combo12 - Combo23     7.70 3.75 6.59   2.052  0.2607
# Combo12 - Combo33     8.60 3.75 6.59   2.292  0.1938
# Combo23 - Combo33     0.90 3.75 6.57   0.240  0.9946
# 
# Treatment = Elevated (18°C):
#   contrast          estimate   SE   df t.ratio p.value
# Combo11 - Combo12     5.40 3.96 8.21   1.362  0.5531
# Combo11 - Combo23    21.20 3.84 7.24   5.521  0.0034
# Combo11 - Combo33    24.82 4.67 7.02   5.319  0.0046
# Combo12 - Combo23    15.80 3.87 7.49   4.080  0.0169
# Combo12 - Combo33    19.42 4.69 7.18   4.138  0.0169
# Combo23 - Combo33     3.62 4.59 6.57   0.788  0.8576
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

tab_model(lmer_pig1)
library(report)
report(lmer_pig1)

library(effects)
e = allEffects(lmer_pig1)
plot(e)

###--------For Dam--------####
lmer_pig_fem=lmer(count ~ Treatment + Dam + (1|Culture:Treatment), data= df)
summary(lmer_pig_fem)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: count ~ Treatment + Dam + (1 | Culture:Treatment)
#    Data: df
# 
# REML criterion at convergence: 2936.9
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.4622 -0.6202 -0.0844  0.5595  4.2162 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 23.03    4.799   
#  Residual                      84.39    9.186   
# Number of obs: 402, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                          Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                55.991      2.242  10.767  24.969 7.07e-11 ***
# TreatmentElevated (18°C)   13.382      2.679  10.998   4.995 0.000406 ***
# Dam2                      -12.490      3.130  10.849  -3.990 0.002180 ** 
# Dam3                      -13.318      3.481  10.748  -3.826 0.002933 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TE(18° Dam2  
# TrtmE(18°C) -0.580              
# Dam2        -0.468 -0.012       
# Dam3        -0.495  0.117  0.305

lmer_pig_fem1=lmer(count ~ Treatment + Dam + Treatment*Dam + (1|Culture:Treatment), data= df)
summary(lmer_pig_fem1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: count ~ Treatment + Dam + Treatment * Dam + (1 | Culture:Treatment)
#    Data: df
# 
# REML criterion at convergence: 2918.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.4718 -0.6310 -0.1045  0.5443  4.1668 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 11.23    3.352   
#  Residual                      84.41    9.188   
# Number of obs: 402, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                               Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                     53.087      1.881   8.296  28.227 1.55e-09 ***
# TreatmentElevated (18°C)        19.409      2.732   9.228   7.105 4.96e-05 ***
# Dam2                            -6.654      3.250   8.218  -2.047   0.0739 .  
# Dam3                            -7.554      3.250   8.218  -2.324   0.0478 *  
# TreatmentElevated (18°C):Dam2  -11.892      4.638   8.522  -2.564   0.0318 *  
# TreatmentElevated (18°C):Dam3  -14.609      5.342   8.436  -2.735   0.0244 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Dam2   Dam3   TE(18°C):D2
# TrtmE(18°C) -0.688                                    
# Dam2        -0.579  0.398                             
# Dam3        -0.579  0.398     0.335                   
# TE(18°C):D2  0.406 -0.589    -0.701 -0.235            
# TE(18°C):D3  0.352 -0.511    -0.204 -0.608  0.301 


anova(lmer_pig_fem1, lmer_pig_fem)
# Data: df
# Models:
#   lmer_pig_fem: count ~ Treatment + Dam + (1 | Culture:Treatment)
# lmer_pig_fem1: count ~ Treatment + Dam + Treatment * Dam + (1 | Culture:Treatment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
# lmer_pig_fem     6 2962.8 2986.7 -1475.4   2950.8                        
# lmer_pig_fem1    8 2954.8 2986.8 -1469.4   2938.8 11.925  2   0.002573 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


em_dam = emmeans(lmer_pig_fem1, specs = pairwise ~ Dam|Treatment)
summary(em_dam)

# Treatment = Ambient (14°C):
#   Dam emmean   SE    df lower.CL upper.CL
# 1        53.1 1.88  8.48     48.8     57.4
# 2        46.4 2.65  8.36     40.4     52.5
# 3        45.5 2.65  8.36     39.5     51.6
# 
# Treatment = Elevated (18°C):
#   Dam emmean   SE    df lower.CL upper.CL
# 1        72.5 1.98 10.43     68.1     76.9
# 2        54.0 2.65  8.36     47.9     60.0
# 3        50.3 3.75  8.36     41.8     58.9
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast          estimate   SE   df t.ratio p.value
# Dam1 - Dam2     6.65 3.25 8.40   2.047  0.1602
# Dam1 - Dam3     7.55 3.25 8.40   2.324  0.1067
# Dam2 - Dam3     0.90 3.75 8.36   0.240  0.9688
# 
# Treatment = Elevated (18°C):
#   contrast          estimate   SE   df t.ratio p.value
# Dam1 - Dam2    18.55 3.31 9.02   5.604  0.0009
# Dam1 - Dam3    22.16 4.24 8.75   5.228  0.0015
# Dam2 - Dam3     3.62 4.59 8.36   0.788  0.7200
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 


##---------Sire-----------##
lmer_pig_sire=lmer(count ~ Treatment + Sire + (1|Culture:Treatment), data= df)
summary(lmer_pig_sire)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: count ~ Treatment + Sire + (1 | Culture:Treatment)
#    Data: df
# 
# REML criterion at convergence: 2936.9
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.4748 -0.6250 -0.0937  0.5596  4.2300 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 22.67    4.762   
#  Residual                      84.39    9.186   
# Number of obs: 402, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                          Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                56.751      2.868  11.026  19.785 5.79e-10 ***
# TreatmentElevated (18°C)   13.441      2.643  11.046   5.085 0.000348 ***
# Sire2                      -1.586      3.629  11.367  -0.437 0.670253    
# Sire3                     -13.631      3.198  11.040  -4.262 0.001327 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TE(18° Sire2 
# TrtmE(18°C) -0.450              
# Sire2       -0.632  0.004       
# Sire3       -0.738  0.049  0.566

lmer_pig_sire1=lmer(count ~ Treatment + Sire + Treatment*Sire + (1|Culture:Treatment), data= df)
summary(lmer_pig_sire1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: count ~ Treatment + Sire + Treatment * Sire + (1 | Culture:Treatment)
#    Data: df
# 
# REML criterion at convergence: 2917.1
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -3.5482 -0.6356 -0.0915  0.5766  4.2632 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept)  9.075   3.012   
#  Residual                      84.383   9.186   
# Number of obs: 402, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                                Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                      52.042      2.454   8.778  21.209 7.61e-09 ***
# TreatmentElevated (18°C)         23.108      3.559   9.715   6.493 7.95e-05 ***
# Sire2                             2.085      3.463   8.702   0.602  0.56246    
# Sire3                            -6.058      2.999   8.704  -2.020  0.07517 .  
# TreatmentElevated (18°C):Sire2   -7.465      5.053   9.863  -1.477  0.17084    
# TreatmentElevated (18°C):Sire3  -16.347      4.427   9.280  -3.692  0.00472 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Sire2  Sire3  TE(18°C):S2
# TrtmE(18°C) -0.689                                    
# Sire2       -0.709  0.489                             
# Sire3       -0.818  0.564     0.580                   
# TE(18°C):S2  0.486 -0.704    -0.685 -0.397            
# TE(18°C):S3  0.554 -0.804    -0.393 -0.677  0.566 


anova(lmer_pig_sire1, lmer_pig_sire)
# 
# refitting model(s) with ML (instead of REML)
# Data: df
# Models:
#   lmer_pig_sire: count ~ Treatment + Sire + (1 | Culture:Treatment)
# lmer_pig_sire1: count ~ Treatment + Sire + Treatment * Sire + (1 | Culture:Treatment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# lmer_pig_sire     6 2962.6 2986.5 -1475.3   2950.6                         
# lmer_pig_sire1    8 2952.4 2984.4 -1468.2   2936.4 14.112  2  0.0008623 ***
  # ---
  # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

em_sire = emmeans(lmer_pig_sire1, specs = pairwise ~ Sire|Treatment)
summary(em_sire)
# 
# $emmeans
# Treatment = Ambient (14°C):
#   Sire emmean   SE    df lower.CL upper.CL
# 1      52.0 2.45  8.72     46.5     57.6
# 2      54.1 2.44  8.57     48.6     59.7
# 3      46.0 1.72  8.51     42.0     49.9
# 
# Treatment = Elevated (18°C):
#   Sire emmean   SE    df lower.CL upper.CL
# 1      75.2 2.58 10.64     69.5     80.8
# 2      69.8 2.63 11.43     64.0     75.5
# 3      52.7 1.99  8.51     48.2     57.3
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast      estimate   SE    df t.ratio p.value
# Sire1 - Sire2    -2.08 3.46  8.65  -0.602  0.8227
# Sire1 - Sire3     6.06 3.00  8.65   2.020  0.1650
# Sire2 - Sire3     8.14 2.99  8.55   2.723  0.0578
# 
# Treatment = Elevated (18°C):
#   contrast      estimate   SE    df t.ratio p.value
# Sire1 - Sire2     5.38 3.68 11.03   1.461  0.3454
# Sire1 - Sire3    22.41 3.26  9.76   6.879  0.0001
# Sire2 - Sire3    17.03 3.30 10.21   5.165  0.0010
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 


#---------------------------------Panelling Figures-----------------------------

#Sire
Sire_pig = ggplot(df, aes(x=Sire, y=count, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black")) + facet_wrap(~Treatment) +
  scale_fill_manual(values=c("#0072B2","darkgrey")) +
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
  #ggtitle("Paternal Effect on Pigment Cell Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Sire") +
  ylab("Pigment Cell Count") +
  scale_y_continuous(limits=c(15,130))

Sire_pig

ggsave("Sire_pig_black_outline.svg", Sire_pig)


#Dam
Dam_pig = ggplot(df, aes(x=Dam, y=count, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black")) + facet_wrap(~Treatment) +
  scale_fill_manual(values=c("#0072B2","darkgrey")) +
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
  #ggtitle("Maternal Effect on Pigment Cell Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Dam") +
  ylab("Pigment Cell Count") +
  scale_y_continuous(limits=c(10,130)) 
Dam_pig

Dam_pig1 = ggplot(df, aes(x=Dam, y=count, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#0072B2","darkgrey")) + 
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
  #ggtitle("Maternal Effect on Pigment Cell Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Dam") +
  ylab("Pigment Cell Count") +
  scale_y_continuous(limits=c(10,130)) 
Dam_pig1



#Full Genotype Figure 
cross_pig = ggplot(df, aes(x=Combo, y=count, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black"))  +
  scale_fill_manual(values=c("#0072B2","darkgrey")) + 
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
  #ggtitle("Maternal Effect on Pigment Cell Count") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Genetic Cross") +
  ylab("Pigment Cell Count") +
  scale_y_continuous(limits=c(10,135)) 

cross_pig


#Environmental Treatment Figure 
treat_pig = ggplot(df, aes(x=Treatment, y=count, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black")) + 
  scale_fill_manual(values=c("#0072B2","darkgrey")) + 
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
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Treatment") +
  ylab("Pigment Cell Count") +
  scale_y_continuous(limits=c(15,120)) +
  geom_signif(comparisons=list(c("Ambient (14°C)", "Elevated (18°C)")), annotations="***", y_position = 111, textsize=7, color="black")


treat_pig


#Genotype Facet
cross_pig = ggplot(df, aes(x=Combo, y=count, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black"))  +
  scale_fill_manual(values=c("#0072B2","darkgrey")) +  facet_wrap(~Treatment) +
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
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Genetic Cross") +
  ylab("Pigment Cell Count") +
  scale_y_continuous(limits=c(15,120))

cross_pig


pigment_panel = ggarrange(treat_pig, cross_pig, Dam_pig, Sire_pig,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2)

pigment_panel


ggsave("pigment_panel_pretty.svg", pigment_panel)

