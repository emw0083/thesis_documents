##Code and stats for correlation between preoral body length and pigment cell count and postoral body length and pigment cell count


#Packages
library(readxl) 
library(readr)
library(tidyr) 
library(reshape2)
library(dplyr)
library(ggplot2) 
theme_set(theme_bw())
cbPalette <- c("#0072B2","darkgrey", "#0072B2","#CC79A7", "#56B4E9", "#E69F00","#999999", "#009E73", "#F0E442", "#D55E00", "#CC79A7")
library(rstatix)
library(lme4)
library(wesanderson)
library(afex)
library(ggpp)
library(ggpmisc)
library(ggpubr)

####--------------- Postoral body length vs pigment cell ------------####
df=read_csv("pigment_and_postoral_for_regression.csv")
#View(df)

#Changes A and E to Ambient (14°C) and Elevated (18°C)
df["Treatment"][df["Treatment"]=="E"]<- "Elevated (18°C)"
df["Treatment"][df["Treatment"]=="A"]<- "Ambient (14°C)"

df_ambient = subset(df, Treatment == "Ambient (14°C)", select =c("culture", "Sample_ID", "stack", "Treatment","tube", "pigment", "avg_mean"))
df_elevated = subset(df, Treatment == "Elevated (18°C)", select =c("culture", "Sample_ID", "stack", "Treatment","tube", "pigment", "avg_mean"))


####--------------Figures---------------####
postoral_all=ggplot(data=df,aes(x=pigment,y=avg_mean, colour=Treatment, fill=Treatment)) +
  geom_smooth(method = "lm", se=TRUE, fill="grey", formula=y~x)+
  labs(y="Pluteus Postoral Body Length (microns)",x="Pluteus Pigment Cell Count")+
  geom_point(size=4, shape=21, color="black") +
  scale_fill_manual(values=(cbPalette)) +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_manual(values=(cbPalette)) + 
  scale_y_continuous(limits=c(150,365)) +
  stat_cor(aes(label=..rr.label..),label.x=95, size=6 )

postoral_all

########-----------------------------------------##########

whole_model=cor.test(df$pigment, df$avg_mean, method = "pearson")
cor.test(df$pigment, df$avg_mean, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  df$pigment and df$avg_mean
# t = 1.3811, df = 399, p-value = 0.168
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.02914775  0.16578775
# sample estimates:
#        cor 
# 0.06897839


report::report(whole_model)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df$pigment and
# df$avg_mean is positive, statistically not significant, and very small
# (r = 0.07, 95% CI [-0.03, 0.17], t(399) = 1.38, p = 0.168)


ambient=cor.test(df_ambient$pigment, df_ambient$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_ambient$pigment and df_ambient$avg_mean
# t = -2.73, df = 233, p-value = 0.006817
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.29732814 -0.04919046
# sample estimates:
#        cor 
# -0.1760546 

report::report(ambient)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df_ambient$pigment and
# df_ambient$avg_mean is negative, statistically significant, and small
# (r = -0.18, 95% CI [-0.30, -0.05], t(233) = -2.73, p = 0.007)


elevated=cor.test(df_elevated$pigment, df_elevated$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_elevated$pigment and df_elevated$avg_mean
# t = -1.9487, df = 164, p-value = 0.05304
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.295974970  0.001930125
# sample estimates:
#        cor 
# -0.1504356

report::report(elevated)

#Effect sizes were labelled following Funder's (2019) recommendations.
#
# The Pearson's product-moment correlation between df_elevated$pigment
# and df_elevated$avg_mean is negative, statistically not significant,
# and small (r = -0.15, 95% CI [-0.30, 1.93e-03], t(164) = -1.95, p =
# 0.053)


ggsave("postoral_regression_updated.svg", postoral_all)



####----------------Preoral Body Length vs pigment cell --------####

df1=read_csv("pigment_preoral_correlation_data.csv")
#View(df1)

#Changes A and E to Ambient (14°C) and Elevated (18°C)
df1["Treatment"][df1["Treatment"]=="E"]<- "Elevated (18°C)"
df1["Treatment"][df1["Treatment"]=="A"]<- "Ambient (14°C)"



####--------------Figures---------------####

preoral_all=ggplot(data=df1,aes(x=pigment,y=avg_mean, colour=Treatment, fill=Treatment)) +
  geom_smooth(method = "lm", se=TRUE, fill="grey", formula=y~x)+
  labs(y="Pluteus Preoral Body Length (microns)",x="Pluteus Pigment Cell Count")+
  geom_point(size=4, shape=21, color="black") +
  scale_fill_manual(values=(cbPalette)) +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_colour_manual(values=(cbPalette)) + stat_cor(aes(label=..rr.label..),label.x=100, size=6 )
 
  
preoral_all


ggsave("preoral_regression_updated.svg", preoral_all)


#####--------------Statistics--------------####

new_stat=lm(avg_mean ~ pigment + Treatment, data=df)
summary(new_stat)

# Call:
#   lm(formula = avg_mean ~ pigment + Treatment, data = df)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -127.368  -14.422   -0.565   14.135   78.290 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)              269.5588     5.2314  51.527  < 2e-16 ***
#   pigment                   -0.3256     0.1009  -3.226  0.00136 ** 
#   TreatmentElevated (18°C)  27.1931     2.7208   9.995  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 23.83 on 398 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.2044,	Adjusted R-squared:  0.2004 
# F-statistic: 51.14 on 2 and 398 DF,  p-value: < 2.2e-16
# 

df$Treatment=as.factor(df$Treatment)

new_stat=lm(avg_mean ~ pigment + relevel(Treatment, ref="Elevated (18°C)"), data=df)
summary(new_stat)




