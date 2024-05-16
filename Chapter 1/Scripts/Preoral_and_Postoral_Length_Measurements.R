####----- Identifying preoral and postoral length in 6dpf larvae reared in two developmental conditions----- ########

#Load Packages
library(readxl)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggplot2)
library("ggsignif")
theme_set(theme_bw())
cbPalette <- c("#0072B2","darkgrey","#84B9EF","#B73E3E", "#DD5353", "#D61C4E", "#FFACAC", "#952E4B", "#0072B2","#CC79A7", "#56B4E9", "#E69F00","#999999", "#009E73", "#F0E442", "#D55E00", "#CC79A7")
library(tidyverse)
library(rstatix)
library(lme4)
library(afex) 
library(ggpubr)
library("gridExtra")
library(dplyr)
library(plyr)
library(multcomp) 
library(emmeans) 

####----------General file manipulation-----------####

raw <- read_excel("Pigment_Cell_and_Length_Data.xlsx")
#View(raw)

#Subset the data so it's easier to work with
data <- raw[ , c("Tube_ID", "Sample_ID", "Photo_ID", "Counter_Type", "LENGTH_MEASUREMENT(micron)" )]
#view(data)

#Subset so it's just by length measurement
types <- subset(data, select = c("Tube_ID", "Sample_ID", "Photo_ID", "Counter_Type", "LENGTH_MEASUREMENT(micron)" ))
#view(types)

#Rename columns so they are easoer to work with
names(types)[names(types) == 'Tube_ID'] <- 'tube'
names(types)[names(types) == 'Sample_ID'] <- 'ID'
names(types)[names(types) == 'Photo_ID'] <- 'stack'
names(types)[names(types) == 'Counter_Type'] <- 'type'
names(types)[names(types) == 'LENGTH_MEASUREMENT(micron)'] <- 'length'

#need to omit rows with NA
df = na.omit(types)
#view(df)

#Aggregates the file so it's just the stack ID, Replicate Culture ID, Culture ID, and mean length value 
means = aggregate(df$length, list(df$stack, df$ID, df$tube), mean)
#view(means)

####------------------------------Preoral Body Length ----------------------####

##Subset data so it's just the preoral values in the file -- there are two preoral types (Preoral1 and Preoral2) since a point was placed on each preoral arm
preoral_mean_sunbset = subset(df, type =="PreOral1" | type == "PreOral2", select =c("tube", "ID", "stack", "type", "length"))

##Takes the average of the PreOral1 and PreOral2 measurements so each individual larvae sampled has one length value
preoral_mean = aggregate(preoral_mean_sunbset$length, list(preoral_mean_sunbset$stack, preoral_mean_sunbset$ID, preoral_mean_sunbset$tube), mean)

##save as a new file for easier accessibility, also is the dataframe for the correlation figure
#write.csv(preoral_mean, "preoral_mean_list.csv")

##Rename columns for identification
colnames(preoral_mean)[colnames(preoral_mean) == "x"] = "mean"
colnames(preoral_mean)[colnames(preoral_mean) == "Group.1"] = "stack"
colnames(preoral_mean)[colnames(preoral_mean) == "Group.2"] = "ID"
colnames(preoral_mean)[colnames(preoral_mean) == "Group.3"] = "tube"

##Extracts values from columns and makes new columns for easier analysis for tested conditions
df$Treatment <- substr(df$tube,2,2)
df$Culture <- substr(df$tube,3,5)
df$Dam <- substr(df$Culture,1,1)
df$Sire <- substr(df$Culture,2,2)
df$Combo <-substr(df$Culture,1,2)

##Extracts values from columns and makes new columns for easier analysis for tested conditions
preoral_mean$Treatment <- substr(preoral_mean$ID,2,2)
preoral_mean$Culture <- substr(preoral_mean$ID,3,5)
preoral_mean$Dam <- substr(preoral_mean$Culture,1,1)
preoral_mean$Sire <- substr(preoral_mean$Culture,2,2)
preoral_mean$Combo <-substr(preoral_mean$Culture,1,2)

##Renames A and E to Ambient(14°C) and Elevated (18°C)
preoral_mean["Treatment"][preoral_mean["Treatment"]=="E"]<- "Elevated (18°C)"
preoral_mean["Treatment"][preoral_mean["Treatment"]=="A"]<- "Ambient (14°C)"


#write.csv(preoral_mean, "preoral_mean_list.csv")

#preoral = read.csv("preoral_mean_list.csv")


####---------Stats for Preoral Body Length ----------------- ####

lmer_all_pre=lmer(mean ~ Treatment + Combo + (1|Culture:Treatment), data = preoral_mean)
summary(lmer_all_pre)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Combo + (1 | Culture:Treatment)
#    Data: preoral_mean
# 
# REML criterion at convergence: 3332.3
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.3802 -0.6315 -0.0233  0.6226  3.1504 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 191.9    13.85   
#  Residual                      230.5    15.18   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                          Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)               249.170      7.995   9.981  31.167 2.81e-11 ***
# TreatmentElevated (18°C)   20.233      7.402   9.983   2.733   0.0211 *  
# Combo12                    -9.764     10.050  10.091  -0.972   0.3540    
# Combo23                    15.168     10.017   9.960   1.514   0.1610    
# Combo33                     3.281     10.882   9.932   0.301   0.7693    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TE(18° Comb12 Comb23
# TrtmE(18°C) -0.459                     
# Combo12     -0.629  0.002              
# Combo23     -0.628 -0.003  0.501       
# Combo33     -0.631  0.111  0.461  0.462


lmer_all_pre1=lmer(mean ~ Treatment + Combo + Treatment*Combo + (1|Culture:Treatment), data = preoral_mean)
summary(lmer_all_pre1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Combo + Treatment * Combo + (1 | Culture:Treatment)
#    Data: preoral_mean
# 
# REML criterion at convergence: 3308.3
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.3704 -0.6370 -0.0288  0.6316  3.1569 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 261.2    16.16   
#  Residual                      230.5    15.18   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                                  Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                       246.155     11.607   6.958  21.207  1.4e-07 ***
# TreatmentElevated (18°C)           26.314     16.465   7.043   1.598    0.154    
# Combo12                            -3.405     16.409   6.947  -0.208    0.842    
# Combo23                            16.071     16.407   6.943   0.980    0.360    
# Combo33                             8.075     16.407   6.943   0.492    0.638    
# TreatmentElevated (18°C):Combo12  -12.879     23.296   7.056  -0.553    0.597    
# TreatmentElevated (18°C):Combo23   -1.856     23.238   6.986  -0.080    0.939    
# TreatmentElevated (18°C):Combo33  -11.417     25.970   6.974  -0.440    0.674    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Comb12 Comb23 Comb33 TE(18°C):C1 TE(18°C):C2
# TrtmE(18°C) -0.705                                                       
# Combo12     -0.707  0.499                                                
# Combo23     -0.707  0.499     0.500                                      
# Combo33     -0.707  0.499     0.500  0.501                               
# TE(18°C):C1  0.498 -0.707    -0.704 -0.353 -0.353                        
# TE(18°C):C2  0.499 -0.709    -0.353 -0.706 -0.353  0.501                 
# TE(18°C):C3  0.447 -0.634    -0.316 -0.316 -0.632  0.448       0.449 

anova(lmer_all_pre1,lmer_all_pre)

# Data: preoral_mean
# Models:
#   lmer_all_pre: mean ~ Treatment + Combo + (1 | Culture:Treatment)
# lmer_all_pre1: mean ~ Treatment + Combo + Treatment * Combo + (1 | Culture:Treatment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# lmer_all_pre     7 3374.4 3402.3 -1680.2   3360.4                     
# lmer_all_pre1   10 3379.5 3419.4 -1679.7   3359.5 0.9144  3     0.8219

eem_all_pre =emmeans(lmer_all_pre1, specs = pairwise ~ Combo|Treatment)
summary(eem_all_pre)

# $emmeans
# Treatment = Ambient (14°C):
#   Combo emmean   SE   df lower.CL upper.CL
# 11       246 11.6 6.96      219      274
# 12       243 11.6 6.94      215      270
# 23       262 11.6 6.93      235      290
# 33       254 11.6 6.93      227      282
# 
# Treatment = Elevated (18°C):
#   Combo emmean   SE   df lower.CL upper.CL
# 11       272 11.7 7.13      245      300
# 12       256 11.7 7.20      229      284
# 23       287 11.6 6.93      259      314
# 33       269 16.4 6.93      230      308
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast          estimate   SE   df t.ratio p.value
# Combo11 - Combo12     3.41 16.4 6.95   0.208  0.9965
# Combo11 - Combo23   -16.07 16.4 6.94  -0.980  0.7653
# Combo11 - Combo33    -8.07 16.4 6.94  -0.492  0.9584
# Combo12 - Combo23   -19.48 16.4 6.93  -1.188  0.6529
# Combo12 - Combo33   -11.48 16.4 6.93  -0.700  0.8939
# Combo23 - Combo33     8.00 16.4 6.93   0.488  0.9594
# 
# Treatment = Elevated (18°C):
#   contrast          estimate   SE   df t.ratio p.value
# Combo11 - Combo12    16.28 16.5 7.17   0.985  0.7624
# Combo11 - Combo23   -14.21 16.5 7.03  -0.864  0.8231
# Combo11 - Combo33     3.34 20.1 6.99   0.166  0.9982
# Combo12 - Combo23   -30.50 16.5 7.07  -1.851  0.3268
# Combo12 - Combo33   -12.94 20.1 7.02  -0.642  0.9150
# Combo23 - Combo33    17.56 20.1 6.93   0.874  0.8181
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates


##-----------------Dam---------------### 
lmer_pre_fem=lmer(mean ~ Treatment + Dam + (1|Culture:Treatment), data= preoral_mean)
summary(lmer_pre_fem)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Dam + (1 | Culture:Treatment)
#    Data: preoral_mean
# 
# REML criterion at convergence: 3339.7
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.3673 -0.6496 -0.0370  0.6225  3.1641 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 190.4    13.80   
#  Residual                      230.5    15.18   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                          Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)               244.288      6.196  10.845  39.430 4.64e-13 ***
# TreatmentElevated (18°C)   20.245      7.375  10.927   2.745   0.0192 *  
# Dam2                       20.044      8.637  10.871   2.321   0.0408 *  
# Dam3                        8.158      9.619  10.834   0.848   0.4147    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TE(18° Dam2  
# TrtmE(18°C) -0.589              
# Dam2        -0.466 -0.004       
# Dam3        -0.494  0.124  0.301

lmer_pre_fem1=lmer(mean ~ Treatment + Dam + Treatment*Dam + (1|Culture:Treatment), data= preoral_mean)
summary(lmer_pre_fem1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Dam + Treatment * Dam + (1 | Culture:Treatment)
#    Data: preoral_mean
# 
# REML criterion at convergence: 3324.1
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.3673 -0.6524 -0.0364  0.6316  3.1615 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 231.1    15.20   
#  Residual                      230.5    15.18   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                               Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                    244.450      7.731   8.875  31.618 1.98e-10 ***
# TreatmentElevated (18°C)        19.905     10.982   9.031   1.813    0.103    
# Dam2                            17.775     13.385   8.858   1.328    0.217    
# Dam3                             9.779     13.385   8.858   0.731    0.484    
# TreatmentElevated (18°C):Dam2    4.554     18.956   8.910   0.240    0.816    
# TreatmentElevated (18°C):Dam3   -5.007     21.880   8.894  -0.229    0.824    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Dam2   Dam3   TE(18°C):D2
# TrtmE(18°C) -0.704                                    
# Dam2        -0.578  0.407                             
# Dam3        -0.578  0.407     0.334                   
# TE(18°C):D2  0.408 -0.579    -0.706 -0.236            
# TE(18°C):D3  0.353 -0.502    -0.204 -0.612  0.291 

anova(lmer_pre_fem1, lmer_pre_fem)
# Data: preoral_mean
# Models:
#   lmer_pre_fem: mean ~ Treatment + Dam + (1 | Culture:Treatment)
# lmer_pre_fem1: mean ~ Treatment + Dam + Treatment * Dam + (1 | Culture:Treatment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# lmer_pre_fem     6 3373.7 3397.7 -1680.9   3361.7                     
# lmer_pre_fem1    8 3377.5 3409.4 -1680.7   3361.5 0.2583  2     0.8788


eem_fem_pre1 =emmeans(lmer_pre_fem1, specs = pairwise ~ Dam|Treatment)
summary(eem_fem_pre1)

# $emmeans
# Treatment = Ambient (14°C):
#   Dam emmean    SE   df lower.CL upper.CL
# 1         244  7.73 8.91      227      262
# 2         262 10.93 8.88      237      287
# 3         254 10.93 8.88      229      279
# 
# Treatment = Elevated (18°C):
#   Dam emmean    SE   df lower.CL upper.CL
# 1         264  7.80 9.22      247      282
# 2         287 10.93 8.88      262      311
# 3         269 15.45 8.88      234      304
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast          estimate   SE   df t.ratio p.value
# Dam1 - Dam2   -17.77 13.4 8.89  -1.328  0.4163
# Dam1 - Dam3    -9.78 13.4 8.89  -0.731  0.7523
# Dam2 - Dam3     8.00 15.5 8.88   0.518  0.8649
# 
# Treatment = Elevated (18°C):
#   contrast          estimate   SE   df t.ratio p.value
# Dam1 - Dam2   -22.33 13.4 8.99  -1.663  0.2705
# Dam1 - Dam3    -4.77 17.3 8.95  -0.276  0.9592
# Dam2 - Dam3    17.56 18.9 8.88   0.928  0.6376
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 


####-----------Sire----------------#### 

lmer_pre_sire=lmer(mean ~ Treatment + Sire + (1|Culture:Treatment), data= preoral_mean)
summary(lmer_pre_sire)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Sire + (1 | Culture:Treatment)
#    Data: preoral_mean
# 
# REML criterion at convergence: 3340.1
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.3787 -0.6426 -0.0228  0.6168  3.1517 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 195.5    13.98   
#  Residual                      230.5    15.18   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                          Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)               248.714      8.056  10.995  30.875 4.92e-12 ***
# TreatmentElevated (18°C)   21.153      7.420  10.994   2.851   0.0158 *  
# Sire2                      -9.764     10.140  11.113  -0.963   0.3561    
# Sire3                      10.135      8.979  10.995   1.129   0.2830    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TE(18° Sire2 
# TrtmE(18°C) -0.457              
# Sire2       -0.629  0.002       
# Sire3       -0.735  0.056  0.564


em_sire1 = emmeans(lmer_pre_sire, specs = pairwise ~ Sire|Treatment)
summary(em_sire1)



lmer_pre_sire1=lmer(mean ~ Treatment + Sire + Sire*Treatment + (1|Culture:Treatment), data= preoral_mean)
summary(lmer_pre_sire1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Sire + Sire * Treatment + (1 | Culture:Treatment)
#    Data: preoral_mean
# 
# REML criterion at convergence: 3324.5
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -4.3711 -0.6380 -0.0279  0.6180  3.1576 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 231.1    15.20   
#  Residual                      230.5    15.18   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                                Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                     246.154     10.940   8.952  22.500 3.46e-09 ***
# TreatmentElevated (18°C)         26.315     15.525   9.075   1.695    0.124    
# Sire2                            -3.405     15.465   8.936  -0.220    0.831    
# Sire3                            12.073     13.394   8.938   0.901    0.391    
# TreatmentElevated (18°C):Sire2  -12.871     21.967   9.094  -0.586    0.572    
# TreatmentElevated (18°C):Sire3   -3.711     19.502   9.014  -0.190    0.853    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Sire2  Sire3  TE(18°C):M2
# TrtmE(18°C) -0.705                                    
# Sire2       -0.707  0.499                             
# Sire3       -0.817  0.576     0.578                   
# TE(18°C):M2  0.498 -0.707    -0.704 -0.407            
# TE(18°C):M3  0.561 -0.796    -0.397 -0.687  0.563     

anova(lmer_pre_sire1, lmer_pre_sire)
# Data: preoral_mean
# Models:
#   lmer_pre_sire: mean ~ Treatment + Sire + (1 | Culture:Treatment)
# lmer_pre_sire1: mean ~ Treatment + Sire + Sire * Treatment + (1 | Culture:Treatment)
# npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# lmer_pre_sire     6 3374.1 3398.0 -1681.0   3362.1                     
# lmer_pre_sire1    8 3377.5 3409.4 -1680.7   3361.5 0.6001  2     0.7408


em_sire_pre1= emmeans(lmer_pre_sire1, specs = pairwise ~ Sire|Treatment)
summary(em_sire_pre1)

# Treatment = Ambient (14°C):
#   Sire emmean    SE   df lower.CL upper.CL
# 1       246 10.94 8.96      221      271
# 2       243 10.93 8.93      218      268
# 3       258  7.73 8.92      241      276
# 
# Treatment = Elevated (18°C):
#   Sire emmean    SE   df lower.CL upper.CL
# 1       272 11.01 9.21      248      297
# 2       256 11.05 9.32      231      281
# 3       281  8.92 8.92      261      301
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast      estimate   SE   df t.ratio p.value
# Sire1 - Sire2     3.40 15.5 8.95   0.220  0.9737
# Sire1 - Sire3   -12.07 13.4 8.95  -0.901  0.6530
# Sire2 - Sire3   -15.48 13.4 8.93  -1.156  0.5064
# 
# Treatment = Elevated (18°C):
#   contrast      estimate   SE   df t.ratio p.value
# Sire1 - Sire2    16.28 15.6 9.26   1.043  0.5695
# Sire1 - Sire3    -8.36 14.2 9.09  -0.590  0.8288
# Sire2 - Sire3   -24.64 14.2 9.16  -1.735  0.2443
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates



####-------------- Figures for Preoral Body Length ----------------####

##By treatment condition 

treat_pre = ggplot(preoral_mean, aes(x=Treatment, y=mean, color=Treatment)) + 
  geom_boxplot() +
  #geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) + 
  geom_point(aes(fill=Treatment, color=Treatment),shape=21, size=2, position = position_jitterdodge(0.25)) +
  scale_color_manual(values=c("black", "black")) + 
  scale_fill_manual(values=c("#0072B2","darkgrey")) + ##this line was previously fill_manual. It needs to match whatever you have in your aes
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
  #ggtitle("Effect of Developmental Temperature on Preoral Body Arm Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(comparisons=list(c("Ambient (14°C)", "Elevated (18°C)")), annotations="***", textsize=7, y_position = 340, color="black") +
  xlab("Treatment") +
  ylab("Preoral Body Length (micron)")

treat_pre

##By Genotype 

genotype_pre =ggplot(preoral_mean, aes(x=Combo, y=mean, color=Treatment)) + 
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
  #ggtitle("Preoral Body Length by Genetic Cross") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Genetic Cross") +
  ylab("Preoral Body Length (micron)") +
  annotate("text", x = 4.40, y = 345, label = "n.s.", size = 6)

genotype_pre


##By Dam 
Dam_pre = ggplot(preoral_mean, aes(x=Dam, y=mean, color=Treatment)) + 
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
  #ggtitle("Maternal Effect on Preoral Body Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Dam") +
  ylab("Preoral Body Length (micron)") +
  annotate("text", x = 3.35, y = 345, label = "n.s.", size = 6)

Dam_pre

##By Sire 

male_pre = ggplot(preoral_mean, aes(x=Sire, y=mean, color=Treatment)) + 
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
  #ggtitle("Paternal Effect on Preoral Body Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Sire") +
  ylab("Preoral Body Length (micron)") +
  annotate("text", x = 3.5, y = 345, label = "n.s.", size = 6)

male_pre



preoral_panel = ggarrange(treat_pre, genotype_pre, Dam_pre, male_pre,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2)

preoral_panel

annotate_figure(preoral_panel,
                top = text_grob("Developmental temperature effects preoral body length", color = "black", face = "bold", size = 14))
                          
ggsave("preoral_length_updated_pretty.svg", preoral_panel)

####--------------------------POSTORAL BODY LENGTH DATA ----------------------####


#For finding postoral body length

##Subset data so it's just the postoral values in the file
postoral_mean_subset = subset(df, type =="PostOral1" | type == "PostOral2", select =c("tube", "ID", "stack", "type", "length"))

##Takes the average of the PostOral1 and PostOral2 measurements so each individual larvae sampled has one length value
postoral_mean = aggregate(postoral_mean_subset$length, list(postoral_mean_subset$stack, postoral_mean_subset$ID, postoral_mean_subset$tube), mean)

##Rename columns so they reflect what they are measuring
colnames(postoral_mean)[colnames(postoral_mean) == "x"] = "mean"
colnames(postoral_mean)[colnames(postoral_mean) == "Group.1"] = "stack"
colnames(postoral_mean)[colnames(postoral_mean) == "Group.2"] = "ID"
colnames(postoral_mean)[colnames(postoral_mean) == "Group.3"] = "tube"

#Makes a file of postoral body length means, this data file is used for the correlation figure 
#write.csv(postoral_mean, "postoral_mean_CORRECT.csv")


##make a column to designate the Elevated vs Ambient
#substr = substring(extracts character values)
##Extracts values from columns and makes new columns for easier analysis for tested conditions
df$Treatment <- substr(df$tube,2,2)
df$Culture <- substr(df$tube,3,5)
df$Dam <- substr(df$Culture,1,1)
df$Sire <- substr(df$Culture,2,2)
df$Combo <-substr(df$Culture,1,2)

##Extracts values from columns and makes new columns for easier analysis for tested conditions
postoral_mean$Treatment <- substr(postoral_mean$ID,2,2)
postoral_mean$Culture <- substr(postoral_mean$ID,3,5)
postoral_mean$Dam <- substr(postoral_mean$Culture,1,1)
postoral_mean$Sire <- substr(postoral_mean$Culture,2,2)
postoral_mean$Combo <-substr(postoral_mean$Culture,1,2)

##Renames A and E to Ambient(14°C) and Elevated (18°C)
postoral_mean["Treatment"][postoral_mean["Treatment"]=="E"]<- "Elevated (18°C)"
postoral_mean["Treatment"][postoral_mean["Treatment"]=="A"]<- "Ambient (14°C)"

####-------------Stats for Postoral Body Length--------------####


#ALL
lmer_all_post=lmer(mean ~ Treatment + Combo + (1|Culture:Treatment), data = postoral_mean)
summary(lmer_all_post)

# Formula: mean ~ Treatment + Combo + (1 | Culture:Treatment)
# Data: postoral_mean
# 
# REML criterion at convergence: 3436.8
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -6.7491 -0.6268  0.0127  0.6531  4.8441 
# 
# Random effects:
#   Groups            Name        Variance Std.Dev.
# Culture:Treatment (Intercept) 282.8    16.82   
# Residual                      299.1    17.30   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#   Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)               247.440      9.679   9.988  25.564 1.96e-10 ***
#   TreatmentElevated (18°C)   22.252      8.962   9.990   2.483   0.0324 *  
#   Combo12                    -7.769     12.163  10.086  -0.639   0.5372    
# Combo23                    11.701     12.128   9.970   0.965   0.3575    
# Combo33                    19.654     13.176   9.945   1.492   0.1668    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#   (Intr) TE(18° Comb12 Comb23
#             TrtmE(18°C) -0.460                     
#             Combo12     -0.628  0.001              
#             Combo23     -0.628 -0.003  0.501       
#             Combo33     -0.630  0.111  0.461  0.462


lmer_all_post1=lmer(mean ~ Treatment + Combo + Treatment*Combo + (1|Culture:Treatment), data = postoral_mean)
summary(lmer_all_post1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Combo + Treatment * Combo + (1 | Culture:Treatment)
#    Data: postoral_mean
# 
# REML criterion at convergence: 3411.7
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -6.7481 -0.6275  0.0143  0.6415  4.8524 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 394.6    19.87   
#  Residual                      299.1    17.30   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                                  Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                       248.376     14.236   6.963  17.447 5.28e-07 ***
# TreatmentElevated (18°C)           20.369     20.185   7.036   1.009    0.346    
# Combo12                            -5.606     20.126   6.953  -0.279    0.789    
# Combo23                             6.925     20.124   6.950   0.344    0.741    
# Combo33                            18.527     20.124   6.950   0.921    0.388    
# TreatmentElevated (18°C):Combo12   -4.407     28.558   7.047  -0.154    0.882    
# TreatmentElevated (18°C):Combo23    9.563     28.496   6.987   0.336    0.747    
# TreatmentElevated (18°C):Combo33    2.456     31.849   6.977   0.077    0.941    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Comb12 Comb23 Comb33 TE(18°C):C1 TE(18°C):C2
# TrtmE(18°C) -0.705                                                       
# Combo12     -0.707  0.499                                                
# Combo23     -0.707  0.499     0.500                                      
# Combo33     -0.707  0.499     0.500  0.500                               
# TE(18°C):C1  0.498 -0.707    -0.705 -0.353 -0.353                        
# TE(18°C):C2  0.500 -0.708    -0.353 -0.706 -0.353  0.501                 
# TE(18°C):C3  0.447 -0.634    -0.316 -0.316 -0.632  0.448       0.449   

em_combo_post = emmeans(lmer_all_post1, specs = pairwise ~ Combo|Treatment)
summary(em_combo_post)

# $emmeans
# Treatment = Ambient (14°C):
#   Combo emmean   SE   df lower.CL upper.CL
# 11       248 14.2 6.96      215      282
# 12       243 14.2 6.94      209      276
# 23       255 14.2 6.94      222      289
# 33       267 14.2 6.94      233      301
# 
# Treatment = Elevated (18°C):
#   Combo emmean   SE   df lower.CL upper.CL
# 11       269 14.3 7.11      235      302
# 12       259 14.3 7.17      225      292
# 23       285 14.2 6.94      252      319
# 33       290 20.1 6.94      242      337
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast          estimate   SE   df t.ratio p.value
# Combo11 - Combo12     5.61 20.1 6.95   0.279  0.9918
# Combo11 - Combo23    -6.93 20.1 6.95  -0.344  0.9848
# Combo11 - Combo33   -18.53 20.1 6.95  -0.921  0.7953
# Combo12 - Combo23   -12.53 20.1 6.94  -0.623  0.9215
# Combo12 - Combo33   -24.13 20.1 6.94  -1.200  0.6462
# Combo23 - Combo33   -11.60 20.1 6.94  -0.577  0.9359
# 
# Treatment = Elevated (18°C):
#   contrast          estimate   SE   df t.ratio p.value
# Combo11 - Combo12    10.01 20.3 7.14   0.494  0.9579
# Combo11 - Combo23   -16.49 20.2 7.02  -0.817  0.8448
# Combo11 - Combo33   -20.98 24.7 7.00  -0.850  0.8296
# Combo12 - Combo23   -26.50 20.2 7.06  -1.312  0.5842
# Combo12 - Combo33   -31.00 24.7 7.02  -1.255  0.6157
# Combo23 - Combo33    -4.50 24.6 6.94  -0.182  0.9976
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

##--------Dam----------##

lmer_fem_post1=lmer(mean ~ Treatment + Dam + Treatment*Dam + (1|Culture:Treatment), data = postoral_mean)
summary(lmer_fem_post1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Dam + Treatment * Dam + (1 | Culture:Treatment)
#    Data: postoral_mean
# 
# REML criterion at convergence: 3427.6
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -6.7448 -0.6177  0.0084  0.6631  4.8429 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 318.8    17.86   
#  Residual                      299.1    17.30   
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                               Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                    245.570      9.072   8.912  27.068 7.26e-10 ***
# TreatmentElevated (18°C)        18.188     12.884   9.060   1.412    0.191    
# Dam2                             9.731     15.707   8.896   0.620    0.551    
# Dam3                            21.333     15.707   8.896   1.358    0.208    
# TreatmentElevated (18°C):Dam2   11.744     22.244   8.945   0.528    0.610    
# TreatmentElevated (18°C):Dam3    4.637     25.674   8.931   0.181    0.861    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Dam2   Dam3   TE(18°C):D2
# TrtmE(18°C) -0.704                                    
# Dam2        -0.578  0.407                             
# Dam3        -0.578  0.407     0.334                   
# TE(18°C):D2  0.408 -0.579    -0.706 -0.236            
# TE(18°C):D3  0.353 -0.502    -0.204 -0.612  0.291   


em_Dam_post = emmeans(lmer_fem_post1, specs = pairwise ~ Dam|Treatment)
summary(em_Dam_post)

# $emmeans
# Treatment = Ambient (14°C):
#   Dam emmean    SE   df lower.CL upper.CL
# 1         246  9.07 8.91      225      266
# 2         255 12.82 8.89      226      284
# 3         267 12.82 8.89      238      296
# 
# Treatment = Elevated (18°C):
#   Dam emmean    SE   df lower.CL upper.CL
# 1         264  9.15 9.21      243      284
# 2         285 12.82 8.89      256      314
# 3         290 18.13 8.89      249      331
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast          estimate   SE   df t.ratio p.value
# Dam1 - Dam2    -9.73 15.7 8.89  -0.620  0.8134
# Dam1 - Dam3   -21.33 15.7 8.89  -1.358  0.4014
# Dam2 - Dam3   -11.60 18.1 8.89  -0.640  0.8026
# 
# Treatment = Elevated (18°C):
#   contrast          estimate   SE   df t.ratio p.value
# Dam1 - Dam2   -21.47 15.8 8.99  -1.363  0.3985
# Dam1 - Dam3   -25.97 20.3 8.95  -1.279  0.4410
# Dam2 - Dam3    -4.50 22.2 8.89  -0.202  0.9777
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 
# 


#####Sire 

lmer_male_post1=lmer(mean ~ Treatment + Sire + Treatment*Sire + (1|Culture:Treatment), data = postoral_mean)
summary(lmer_male_post1)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: mean ~ Treatment + Sire + Treatment * Sire + (1 | Culture:Treatment)
#    Data: postoral_mean
# 
# REML criterion at convergence: 3428
# 
# Scaled residuals: 
#     Min      1Q  Median      3Q     Max 
# -6.7498 -0.6259  0.0164  0.6499  4.8404 
# 
# Random effects:
#  Groups            Name        Variance Std.Dev.
#  Culture:Treatment (Intercept) 320.5    17.9    
#  Residual                      299.1    17.3    
# Number of obs: 401, groups:  Culture:Treatment, 15
# 
# Fixed effects:
#                                Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)                     248.375     12.869   8.947  19.300 1.34e-08 ***
# TreatmentElevated (18°C)         20.370     18.258   9.062   1.116    0.293    
# Sire2                            -5.604     18.192   8.932  -0.308    0.765    
# Sire3                            12.727     15.756   8.934   0.808    0.440    
# TreatmentElevated (18°C):Sire2   -4.397     25.834   9.080  -0.170    0.869    
# TreatmentElevated (18°C):Sire3    5.259     22.938   9.005   0.229    0.824    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#             (Intr) TrE(18°C) Sire2  Sire3  TE(18°C):M2
# TrtmE(18°C) -0.705                                    
# Sire2       -0.707  0.499                             
# Sire3       -0.817  0.576     0.578                   
# TE(18°C):M2  0.498 -0.707    -0.704 -0.407            
# TE(18°C):M3  0.561 -0.796    -0.397 -0.687  0.563

em_male_post = emmeans(lmer_male_post1, specs = pairwise ~ Sire|Treatment)
summary(em_male_post)

# emmeans
# Treatment = Ambient (14°C):
#   Sire emmean    SE   df lower.CL upper.CL
# 1       248 12.87 8.96      219      278
# 2       243 12.86 8.93      214      272
# 3       261  9.09 8.92      241      282
# 
# Treatment = Elevated (18°C):
#   Sire emmean    SE   df lower.CL upper.CL
# 1       269 12.95 9.20      240      298
# 2       259 12.99 9.30      230      288
# 3       287 10.50 8.92      263      311
# 
# Degrees-of-freedom method: kenward-roger 
# Confidence level used: 0.95 
# 
# $contrasts
# Treatment = Ambient (14°C):
#   contrast      estimate   SE   df t.ratio p.value
# Sire1 - Sire2      5.6 18.2 8.95   0.308  0.9493
# Sire1 - Sire3    -12.7 15.8 8.95  -0.808  0.7079
# Sire2 - Sire3    -18.3 15.7 8.93  -1.164  0.5021
# 
# Treatment = Elevated (18°C):
#   contrast      estimate   SE   df t.ratio p.value
# Sire1 - Sire2     10.0 18.3 9.25   0.545  0.8513
# Sire1 - Sire3    -18.0 16.7 9.09  -1.079  0.5494
# Sire2 - Sire3    -28.0 16.7 9.15  -1.676  0.2651
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 3 estimates 



####---------------- Figures for Postoral Body length --------------####

##By Environment
treat_post = ggplot(postoral_mean, aes(x=Treatment, y=mean, color=Treatment)) + 
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
  #ggtitle("Effect of Developmental Temperature on Postoral Body Arm Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Treatment") +
  ylab("Postoral Body Length (micron)") +
  geom_signif(comparisons=list(c("Ambient (14°C)", "Elevated (18°C)")), annotations="***", textsize=7, y_position = 352, color="black") +
  scale_y_continuous(limits=c(165,365))
treat_post


##By genotype
cross_post = ggplot(postoral_mean, aes(x=Combo, y=mean, color=Treatment)) + 
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
  #ggtitle("Postoral Body Length by Genetic Cross") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Genetic Cross") +
  ylab("Postoral Body Length (micron)") +
  annotate("text", x = 4.4, y = 365, label = "n.s.", size = 7) +
  scale_y_continuous(limits=c(165,365))

cross_post


##By Dam 
Dam_post = ggplot(postoral_mean, aes(x=Dam, y=mean, color=Treatment)) + 
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
  #ggtitle("Maternal Effect on Postoral Body Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Dam") +
  ylab("Postoral Body Length (micron)") +
  scale_y_continuous(limits=c(165,365)) +
  annotate("text", x = 3.5, y = 365, label = "n.s.", size = 7)

Dam_post

##By Sire 
male_post = ggplot(postoral_mean, aes(x=Sire, y=mean, color=Treatment)) + 
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
  #ggtitle("Paternal Effect on Postoral Body Length") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Sire") +
  ylab("Postoral Body Length (micron)") +
  scale_y_continuous(limits=c(165,365)) + 
  annotate("text", x = 3.5, y = 365, label = "n.s.", size = 7)

male_post



##Panels All Figures 
postoral_panel = ggarrange(treat_post, cross_post, Dam_post, male_post,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2)

postoral_panel

#Adds title to paneled figure
annotate_figure(postoral_panel,
                top = text_grob("Developmental temperature effects postoral body length", color = "black", face = "bold", size = 14))


ggsave("postoral_length_updated_pretty.svg", postoral_panel)


