library(ggplot2)
library(tidyverse)

setwd("~/Desktop/Phadnis_Lab/human_inversions/hwe_scan/are_SDs_out_of_HWE")


#openning files
win_to_NoSD_df <- read.table("windows_NOTintersecting_SDs.bed", header=FALSE, col.names=c("CHROM", "START", "END", "P_PROP"))
win_to_SD_df <- read.table("windows_intersecting_SDs.bed", header=FALSE, col.names=c("CHROM", "START", "END", "P_PROP"))
win_to_LargeSD_df <- read.table("windows_intersecting_LargeSDs.bed", header=FALSE, col.names=c("CHROM", "START", "END", "P_PROP"))

win_mapped_to_SDs_df <- read.table("windows_mapped_to_SDs.bed", header=FALSE, col.names=c("CHROM", "START", "END", "P_PROP", "SDCHROM", "SDSTART", "SDEND"))
win_overlap_with_SDs_df <- read.table("windows_SD_overlaps.bed", header=FALSE, col.names=c("CHROM", "START", "END", "P_PROP", "SDCHROM", "SDSTART", "SDEND", "OVERLAP"))

#label features
win_to_NoSD_df$FEATURE <- rep("No SDs", nrow(win_to_NoSD_df))
win_to_SD_df$FEATURE <- rep("All SDs", nrow(win_to_SD_df))
win_to_LargeSD_df$FEATURE <- rep("Large SDs", nrow(win_to_LargeSD_df))

#combine intersect dfs
windows_SDs_PProps <- rbind(win_to_NoSD_df, win_to_SD_df, win_to_LargeSD_df)
windows_SDs_PProps$FEATURE <- factor(windows_SDs_PProps$FEATURE, levels=c("No SDs", "All SDs", "Large SDs"))



############################################################################################


#anova
linear_model <- lm(P_PROP ~ FEATURE, data=windows_SDs_PProps)
summary(linear_model)
# Residuals:
#     Min       1Q   Median       3Q      Max 
# -0.20185 -0.02151 -0.00296  0.01446  0.92555 
# 
# Coefficients:
#     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.0744475  0.0001172   635.1   <2e-16 ***
#     FEATUREAll SDs   0.0334258  0.0002388   140.0   <2e-16 ***
#     FEATURELarge SDs 0.1273985  0.0005550   229.5   <2e-16 ***
# 
# Residual standard error: 0.05525 on 303006 degrees of freedom
# Multiple R-squared:  0.1796,	Adjusted R-squared:  0.1796 
# F-statistic: 3.316e+04 on 2 and 303006 DF,  p-value: < 2.2e-16
model_anova <- aov(linear_model)
summary(model_anova)
# Df Sum Sq Mean Sq F value Pr(>F)    
# FEATURE          2  202.5   101.2   33165 <2e-16 ***
#     Residuals   303006  925.0     0.0                   
tukeytest <- TukeyHSD(model_anova)
tukeytest
# Tukey multiple comparisons of means
# 95% family-wise confidence level
# 
# Fit: aov(formula = linear_model)
# 
# $FEATURE
# diff        lwr        upr p adj
# All SDs-No SDs    0.03342582 0.03286606 0.03398559     0
# Large SDs-No SDs  0.12739847 0.12609771 0.12869924     0
# Large SDs-All SDs 0.09397265 0.09261090 0.09533440     0


res = 300
scale = 1.5
png("P_density_boxplot.png", height = 4 * scale * res, width = 3 * scale * res, res = res)
p <- ggplot(data=windows_SDs_PProps, aes(x=FEATURE, y=P_PROP, fill=FEATURE)) +
    geom_boxplot(outliers=FALSE, show.legend = FALSE) +
    scale_fill_manual(values=c("#66C2A5", "#8DA0CB", "#FC8D62")) +
    ylab("Proportion of HWE p < 0.01") + xlab("") +
    theme_bw()
print(p)
dev.off()


res = 300
scale = 2
png("Tukey_post_hoc.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
p <- plot(tukeytest)
print(p)
dev.off()


############################################################################################


#filter for windows that don't map to SDs
win_mapped_to_SDs_df <- win_mapped_to_SDs_df[win_mapped_to_SDs_df$SDCHROM != ".", ]

#calculate SD lengths
win_mapped_to_SDs_df$SDLENGTH <- win_mapped_to_SDs_df$SDEND - win_mapped_to_SDs_df$SDSTART


#linear regression
pprop_vs_SDlen <- lm(P_PROP ~ SDLENGTH, data=win_mapped_to_SDs_df)
summary(pprop_vs_SDlen)
# Residuals:
#     Min       1Q   Median       3Q      Max 
# -0.32758 -0.03321 -0.01576  0.00744  0.89509 
# 
# Coefficients:
#     Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 9.786e-02  3.084e-04   317.4   <2e-16 ***
#     SDLENGTH    1.345e-07  1.203e-09   111.8   <2e-16 ***
# 
# Residual standard error: 0.07835 on 70496 degrees of freedom
# Multiple R-squared:  0.1505,	Adjusted R-squared:  0.1505 
# F-statistic: 1.249e+04 on 1 and 70496 DF,  p-value: < 2.2e-16
anova_on_model <- anova(pprop_vs_SDlen)
anova_df <- data.frame(anova_on_model)
anova_df <- anova_df %>%
    select(Sum.Sq)
summedsumsq <- sum(anova_df$Sum.Sq)
anova_df$PVE <- anova_df$Sum.Sq / summedsumsq
print(anova_df)
# Sum.Sq       PVE
# SDLENGTH   76.67776 0.1505279
# Residuals 432.71464 0.8494721


res = 300
scale = 1.5
png("Pden_vs_SDlen.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
p <- ggplot(data=win_mapped_to_SDs_df, aes(x=SDLENGTH, y=P_PROP)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method="lm") +
    ylab("Proportion of HWE p < 0.01") + xlab("Length (bp) of SD that intersects window") +
    theme_bw()
print(p)
dev.off()



############################################################################################

#filter for windows that don't map to SDs
win_overlap_with_SDs_df <- win_overlap_with_SDs_df[win_overlap_with_SDs_df$SDCHROM != ".", ]


#linear regression
pprop_vs_SDoverlap <- lm(P_PROP ~ OVERLAP, data=win_overlap_with_SDs_df)
summary(pprop_vs_SDoverlap)
# Residuals:
#     Min       1Q   Median       3Q      Max 
# -0.22769 -0.02827 -0.00866  0.01301  0.91578 
# 
# Coefficients:
#     Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 8.324e-02  3.225e-04   258.1   <2e-16 ***
#     OVERLAP     1.445e-06  9.559e-09   151.1   <2e-16 ***
# 
# Residual standard error: 0.07388 on 70496 degrees of freedom
# Multiple R-squared:  0.2447,	Adjusted R-squared:  0.2447 
# F-statistic: 2.284e+04 on 1 and 70496 DF,  p-value: < 2.2e-16
anova_on_model <- anova(pprop_vs_SDoverlap)
anova_df <- data.frame(anova_on_model)
anova_df <- anova_df %>%
    select(Sum.Sq)
summedsumsq <- sum(anova_df$Sum.Sq)
anova_df$PVE <- anova_df$Sum.Sq / summedsumsq
print(anova_df)
# Sum.Sq       PVE
# OVERLAP   124.6365 0.2446769
# Residuals 384.7559 0.7553231



res = 300
scale = 1.5
png("Pden_vs_SDoverlap.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
p <- ggplot(data=win_overlap_with_SDs_df, aes(x=OVERLAP, y=P_PROP)) +
    geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method="lm") +
    ylab("Proportion of HWE p < 0.01") + xlab("Length (bp) of Overlap between SD and window") +
    theme_bw()
print(p)
dev.off()




