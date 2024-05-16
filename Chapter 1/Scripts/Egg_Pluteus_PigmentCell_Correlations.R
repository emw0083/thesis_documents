###correlation between egg diameter and larval length

###going to have egg diameter on x-axis and length/pigment cell on y

#Packages
#Packages
library(readxl) #lets us import an excel file
library(readr)
library(tidyr) 
library(reshape2)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
cbPalette <- c("#643d6e", "#a982b4", "#d4c0d9")

library(rstatix)
library(lme4)
library(wesanderson)
library(afex)
library(ggpp)
library(ggpmisc)
library(ggpubr)
library(MuMIn)

####--------------
#read in data file

df=read_csv("egg_postoral_ambient.csv")
#View(df)

#Changes A and E to Ambient (14°C) and Elevated (18°C)
df["Treatment"][df["Treatment"]=="A"]<- "Ambient (14°C)"
df = na.omit(df)

##subset by dam 
df_fem1 = subset(df, Female == "1", select =c("ID", "tube...2", "avg_diameter", "culture","Treatment", "tube...6", "pigment", "avg_mean", "Female"))
df_fem2 = subset(df, Female == "2", select =c("ID", "tube...2", "avg_diameter", "culture","Treatment", "tube...6", "pigment", "avg_mean", "Female"))
df_fem3 = subset(df, Female == "3", select =c("ID", "tube...2", "avg_diameter", "culture","Treatment", "tube...6", "pigment", "avg_mean", "Female"))

df$Female=as.factor(df$Female)
#df$avg_diameter=as.factor(df$avg_diameter)
#df$avg_mean=as.factor(df$avg_mean)


ambient1=ggplot(data=df,aes(x=avg_diameter,y=avg_mean, colour=Female, fill=Female)) +
  geom_smooth(method = "lm", se=TRUE, fill="grey", formula=y~x)+
  labs(y="Pluteus Postoral Body Length (microns)",x="Egg Diameter")+
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
  stat_cor(aes(label=..rr.label..),label.x=98, size=4 )
  #stat_cor(aes(label=..p.value..),label.x=80, size=4 )
  # stat_cor(
  #   aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
  #   label.x = 98
  # )
  
ambient1

ggsave("postoral_egg_regression_ambient.svg", ambient1)


#######---------Stats------------------######

whole_model=cor.test(df$avg_diameter, df$avg_mean, method = "pearson")
cor.test(df$avg_diameter, df$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df$avg_diameter and df$avg_mean
# t = -2.9752, df = 178, p-value = 0.003334
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.35269236 -0.07374134
# sample estimates:
#        cor 
# -0.2176571 

report::report(whole_model)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df$avg_diameter and
# df$avg_mean is negative, statistically significant, and medium 
#(r =-0.22, 95% CI [-0.35, -0.07], t(178) = -2.98, p = 0.003)
# 


female1=cor.test(df_fem1$avg_diameter, df_fem1$avg_mean, method = "pearson")
cor.test(df_fem1$avg_diameter, df_fem1$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# # data:  df_fem1$avg_diameter and df_fem1$avg_mean
# # t = 0.48591, df = 58, p-value = 0.6289
# # alternative hypothesis: true correlation is not equal to 0
# # 95 percent confidence interval:
# #  -0.1933775  0.3125452
# # sample estimates:
# #        cor 
# # 0.06367382

report::report(female1)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df_fem1$avg_diameter
# and df_fem1$avg_mean is positive, statistically not significant, and
# very small (r = 0.06, 95% CI [-0.19, 0.31], t(58) = 0.49, p = 0.629)
# 
# 


female2=cor.test(df_fem2$avg_diameter, df_fem2$avg_mean, method = "pearson")
cor.test(df_fem2$avg_diameter, df_fem2$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_fem2$avg_diameter and df_fem2$avg_mean
# t = 0.40909, df = 58, p-value = 0.684
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2030510  0.3034312
# sample estimates:
#        cor 
# 0.05363931 


report::report(female2)

# The Pearson's product-moment correlation between df_fem2$avg_diameter
# and df_fem2$avg_mean is positive, statistically not significant, and
# very small (r = 0.05, 95% CI [-0.20, 0.30], t(58) = 0.41, p = 0.684)
# 



female3=cor.test(df_fem3$avg_diameter, df_fem3$avg_mean, method = "pearson")
cor.test(df_fem3$avg_diameter, df_fem3$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_fem3$avg_diameter and df_fem3$avg_mean
# t = 0.34145, df = 58, p-value = 0.734
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2115409  0.2953553
# sample estimates:
#       cor 
# 0.0447897 

report::report(female3)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df_fem3$avg_diameter
# and df_fem3$avg_mean is positive, statistically not significant, and
# tiny (r = 0.04, 95% CI [-0.21, 0.30], t(58) = 0.34, p = 0.734)


ggsave("postoral_egg_regression_ambient.svg", ambient1)


####---------------pigment cell ambient------------####

ambientp=ggplot(data=df,aes(x=avg_diameter,y=pigment, colour=Female, fill=Female)) +
  geom_smooth(method = "lm", se=TRUE, fill="grey", formula=y~x)+
  labs(y="Pigment Cell Count",x="Egg Diameter")+
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
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 98
  )

ambientp

ggsave("pigment_egg_regression_ambient.svg", ambientp)



###-----------------Statistics-------------------####
whole_model=cor.test(df$avg_diameter, df$pigment, method = "pearson")
cor.test(df$avg_diameter, df$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df$avg_diameter and df$pigment
# t = 3.635, df = 178, p-value = 0.0003636
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1212747 0.3939899
# sample estimates:
#       cor 
# 0.2628751 

report::report(whole_model)

# The Pearson's product-moment correlation between df$avg_diameter and
# df$pigment is positive, statistically significant, and medium (r =
# 0.26, 95% CI [0.12, 0.39], t(178) = 3.64, p < .001)


female1=cor.test(df_fem1$avg_diameter, df_fem1$pigment, method = "pearson")
cor.test(df_fem1$avg_diameter, df_fem1$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_fem1$avg_diameter and df_fem1$pigment
# t = 0.24024, df = 58, p-value = 0.811
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2241898  0.2831872
# sample estimates:
#        cor 
# 0.03152977

report::report(female1)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df_fem1$avg_diameter
# and df_fem1$pigment is positive, statistically not significant, and
# tiny (r = 0.03, 95% CI [-0.22, 0.28], t(58) = 0.24, p = 0.811)



female2=cor.test(df_fem2$avg_diameter, df_fem2$pigment, method = "pearson")
cor.test(df_fem2$avg_diameter, df_fem2$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_fem2$avg_diameter and df_fem2$pigment
# t = -0.32582, df = 58, p-value = 0.7457
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.2934820  0.2134993
# sample estimates:
#         cor 
# -0.04274262

report::report(female2)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df_fem2$avg_diameter
# and df_fem2$pigment is negative, statistically not significant, and
# tiny (r = -0.04, 95% CI [-0.29, 0.21], t(58) = -0.33, p = 0.746)


female3=cor.test(df_fem3$avg_diameter, df_fem3$pigment, method = "pearson")
cor.test(df_fem3$avg_diameter, df_fem3$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df_fem3$avg_diameter and df_fem3$pigment
# t = 2.2964, df = 58, p-value = 0.02529
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.03751693 0.50555590
# sample estimates:
#       cor 
# 0.2886914 

report::report(female3)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df_fem3$avg_diameter
# and df_fem3$pigment is positive, statistically significant, and medium
# (r = 0.29, 95% CI [0.04, 0.51], t(58) = 2.30, p = 0.025)




########-------------elevated--------------#######

df1=read_csv("egg_postoral_elevated.csv")
#View(df)

#Changes A and E to Ambient (14°C) and Elevated (18°C)
df1["Treatment"][df1["Treatment"]=="E"]<- "Elevated (18°C)"

df1$Female=as.factor(df1$Female)
df1 = na.omit(df1)

##subset by female 
df1_fem1 = subset(df1, Female == "1", select =c("ID", "tube...2", "avg_diameter", "culture","Treatment", "tube...9", "pigment", "avg_mean", "Female"))
df1_fem2 = subset(df1, Female == "2", select =c("ID", "tube...2", "avg_diameter", "culture","Treatment", "tube...9", "pigment", "avg_mean", "Female"))
df1_fem3 = subset(df1, Female == "3", select =c("ID", "tube...2", "avg_diameter", "culture","Treatment", "tube...9", "pigment", "avg_mean", "Female"))

elevated1=ggplot(data=df1,aes(x=avg_diameter,y=avg_mean, colour=Female, fill=Female)) +
  geom_smooth(method = "lm", se=TRUE, fill="grey", formula=y~x)+
  labs(y="Pluteus Postoral Body Length (microns)",x="Egg Diameter")+
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
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 96
  )

elevated1


ggsave("postoral_egg_regression_elevated.svg", elevated1)

########---------------STATS--------------###########
whole_model=cor.test(df1$avg_diameter, df1$avg_mean, method = "pearson")
cor.test(df1$avg_diameter, df1$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1$avg_diameter and df1$avg_mean
# t = -1.0305, df = 88, p-value = 0.3056
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3092956  0.1001572
# sample estimates:
#        cor 
# -0.1091989

report::report(whole_model)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1$avg_diameter and
# df1$avg_mean is negative, statistically not significant, and small (r =
# -0.11, 95% CI [-0.31, 0.10], t(88) = -1.03, p = 0.306)


female1=cor.test(df1_fem1$avg_diameter, df1_fem1$avg_mean, method = "pearson")
cor.test(df1_fem1$avg_diameter, df1_fem1$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1_fem1$avg_diameter and df1_fem1$avg_mean
# t = 0.98541, df = 28, p-value = 0.3329
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1897041  0.5097265
# sample estimates:
#       cor 
# 0.1830774 

report::report(female1)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1_fem1$avg_diameter
# and df1_fem1$avg_mean is positive, statistically not significant, and
# small (r = 0.18, 95% CI [-0.19, 0.51], t(28) = 0.99, p = 0.333)
# 
# 


female2=cor.test(df1_fem2$avg_diameter, df1_fem2$avg_mean, method = "pearson")
cor.test(df1_fem2$avg_diameter, df1_fem2$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1_fem2$avg_diameter and df1_fem2$avg_mean
# t = -0.013335, df = 28, p-value = 0.9895
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3624603  0.3580742
# sample estimates:
#          cor 
# -0.002520124

report::report(female2)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1_fem2$avg_diameter
# and df1_fem2$avg_mean is negative, statistically not significant, and
# tiny (r = -2.52e-03, 95% CI [-0.36, 0.36], t(28) = -0.01, p = 0.989)
# 

female3=cor.test(df1_fem3$avg_diameter, df1_fem3$avg_mean, method = "pearson")
cor.test(df1_fem3$avg_diameter, df1_fem3$avg_mean, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1_fem3$avg_diameter and df1_fem3$avg_mean
# t = 2.7043, df = 28, p-value = 0.01151
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1133940 0.7004963
# sample estimates:
#       cor 
# 0.4550724 

report::report(female3)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1_fem3$avg_diameter
# and df1_fem3$avg_mean is positive, statistically significant, and very
# large (r = 0.46, 95% CI [0.11, 0.70], t(28) = 2.70, p = 0.012)
# 



####------------------pigment cell for elevated-------------####

elevatedp=ggplot(data=df1,aes(x=avg_diameter,y=pigment, colour=Female, fill=Female)) +
  geom_smooth(method = "lm", se=TRUE, fill="grey", formula=y~x)+
  labs(y="Pigment Cell Count",x="Egg Diameter")+
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
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    label.x = 96
  )

elevatedp

ggsave("pigment_egg_regression_elevated.svg", elevatedp)

#####----------------Stats-----------------#####
whole_model=cor.test(df1$avg_diameter, df1$pigment, method = "pearson")
cor.test(df1$avg_diameter, df1$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1$avg_diameter and df1$pigment
# t = 3.8107, df = 88, p-value = 0.0002563
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1835697 0.5412591
# sample estimates:
#       cor 
# 0.3763535

report::report(whole_model)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1$avg_diameter and df1$pigment
# is positive, statistically significant, and large (r = 0.38, 95% CI [0.18, 0.54],
# t(88) = 3.81, p < .001)
# 


female1=cor.test(df1_fem1$avg_diameter, df1_fem1$pigment, method = "pearson")
cor.test(df1_fem1$avg_diameter, df1_fem1$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1_fem1$avg_diameter and df1_fem1$pigment
# t = -0.28336, df = 28, p-value = 0.779
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.4059223  0.3128226
# sample estimates:
#         cor 
# -0.05347306 

report::report(female1)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1_fem1$avg_diameter and
# df1_fem1$pigment is negative, statistically not significant, and very small 
#(r =-0.05, 95% CI [-0.41, 0.31], t(28) = -0.28, p = 0.779)

female2=cor.test(df1_fem2$avg_diameter, df1_fem2$pigment, method = "pearson")
cor.test(df1_fem2$avg_diameter, df1_fem2$pigment, method = "pearson")

# Pearson's product-moment correlation
# 
# data:  df1_fem2$avg_diameter and df1_fem2$pigment
# t = 0.35557, df = 28, p-value = 0.7248
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.3004815  0.4172368
# sample estimates:
#        cor 
# 0.06704575

report::report(female2)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1_fem2$avg_diameter and
# df1_fem2$pigment is positive, statistically not significant, and very small 
#(r = 0.07, 95% CI [-0.30, 0.42], t(28) = 0.36, p = 0.725)


female3=cor.test(df1_fem3$avg_diameter, df1_fem3$pigment, method = "pearson")
cor.test(df1_fem3$avg_diameter, df1_fem3$pigment, method = "pearson")


# Pearson's product-moment correlation
# 
# data:  df1_fem3$avg_diameter and df1_fem3$pigment
# t = 0.010404, df = 28, p-value = 0.9918
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.358557  0.361979
# sample estimates:
#         cor 
# 0.001966209 


report::report(female3)

# Effect sizes were labelled following Funder's (2019) recommendations.
# 
# The Pearson's product-moment correlation between df1_fem3$avg_diameter and
# df1_fem3$pigment is positive, statistically not significant, and tiny 
#(r =1.97e-03, 95% CI [-0.36, 0.36], t(28) = 0.01, p = 0.992)

