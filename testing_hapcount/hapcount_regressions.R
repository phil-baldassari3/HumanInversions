#version 3 of this script that regresses unique haplotype proportion against breakpoint vs. non-breakpoint and SD non-SD

library(ggplot2)
library(data.table)
library(car)
library(tidyverse)



categorical_boxplot <- function(xdata, ydata, plottitle, xlabel, ylabel, outfile) {
  
  # #linear model
  # model <- lm(ydata ~ xdata)
  # anova_on_model <- anova(model)
  # anova_df <- data.frame(anova_on_model)
  # anova_df <- anova_df %>%
  #   select(Sum.Sq)
  # summedsumsq <- sum(anova_df$Sum.Sq)
  # anova_df$PVE <- anova_df$Sum.Sq / summedsumsq
  # print(anova_df)
  
  #making ploting df
  df <- data.frame(xcol = xdata, ycol = ydata)
  
  #wilcoxon test
  categories <- unique(df$xcol)
  dfcat1 <- df[df$xcol == categories[1], ]
  dfcat2 <- df[df$xcol == categories[2], ]
  result <- wilcox.test(dfcat1$ycol, dfcat2$ycol)
  print(result)
  
  #plotting
  res = 300
  scale = 1.3
  png(outfile, height = 4 * scale * res, width = 3 * scale * res, res = res)
  p <- ggplot(df, aes(x=xcol, y=ycol, fill=xcol)) +
    geom_boxplot(outliers = FALSE, width = 0.5) +
    #geom_smooth(method='lm', aes(group = 1), se=TRUE) + 
    labs(x=xlabel, y=ylabel, title=plottitle) +
    theme_classic() +
    theme(legend.position = "none", plot.title = element_text(size = 10))
  print(p)
  dev.off()
  
  res = 300
  scale = 1.3
  png(paste("density_", outfile), height = 3 * scale * res, width = 4 * scale * res, res = res)
  p2 <- ggplot(df, aes(x=ycol, fill=xcol)) +
    geom_density(position = "identity", alpha = 0.5, linewidth = 0.1) +
    scale_fill_discrete(name = "Window type") +
    labs(x=ylabel, title=plottitle) +
    theme_classic()
  print(p2)
  dev.off()
  
  
}



setwd("~/Desktop/Phadnis_Lab/human_inversions/testing_hapcount_classification_power")


d1 <- "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount_labeled.tsv"
d2 <- "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount_labeled.tsv"
d3 <- "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount_labeled.tsv"
d4 <- "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount_labeled.tsv"



###10Kb###
df <- as.data.frame(fread(d1, sep="\t", header=TRUE))
df$InvBrk <- as.factor(df$InvBrk)
df$LargeInvBrk <- as.factor(df$LargeInvBrk)
df$SD <- factor(df$SD, level=c("SD", "Non-SD"))
df$LargeSD <- factor(df$LargeSD, level=c("LargeSD", "Non-LargeSD"))

###INV BRK###
categorical_boxplot(df$InvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nBreakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.brkpnt_10Kb.png")
###LARGE INV BRK###
categorical_boxplot(df$LargeInvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge Inv Breakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.largebrkpnt_10Kb.png")
###SD###
categorical_boxplot(df$SD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nSD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.SD_10Kb.png")
###Large SD###
categorical_boxplot(df$LargeSD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge SD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.largeSD_10Kb.png")



###############################################################
###############################################################


###100Kb###
df <- as.data.frame(fread(d2, sep="\t", header=TRUE))
df$InvBrk <- as.factor(df$InvBrk)
df$LargeInvBrk <- as.factor(df$LargeInvBrk)
df$SD <- factor(df$SD, level=c("SD", "Non-SD"))
df$LargeSD <- factor(df$LargeSD, level=c("LargeSD", "Non-LargeSD"))


###INV BRK###
categorical_boxplot(df$InvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nBreakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.brkpnt_100Kb.png")
###LARGE INV BRK###
categorical_boxplot(df$LargeInvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge Inv Breakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.largebrkpnt_100Kb.png")
###SD###
categorical_boxplot(df$SD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nSD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.SD_100Kb.png")
###Large SD###
categorical_boxplot(df$LargeSD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge SD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.largeSD_100Kb.png")


###############################################################
###############################################################


###1KSNP###
df <- as.data.frame(fread(d3, sep="\t", header=TRUE))
df$InvBrk <- as.factor(df$InvBrk)
df$LargeInvBrk <- as.factor(df$LargeInvBrk)
df$SD <- factor(df$SD, level=c("SD", "Non-SD"))
df$LargeSD <- factor(df$LargeSD, level=c("LargeSD", "Non-LargeSD"))


###INV BRK###
categorical_boxplot(df$InvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nBreakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.brkpnt_1KSNP.png")
###LARGE INV BRK###
categorical_boxplot(df$LargeInvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge Inv Breakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.largebrkpnt_1KSNP.png")
###SD###
categorical_boxplot(df$SD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nSD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.SD_10Kb.png")
###Large SD###
categorical_boxplot(df$LargeSD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge SD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.largeSD_1KSNP.png")


###############################################################
###############################################################


###10KSNP###
df <- as.data.frame(fread(d4, sep="\t", header=TRUE))
df$InvBrk <- as.factor(df$InvBrk)
df$LargeInvBrk <- as.factor(df$LargeInvBrk)
df$SD <- factor(df$SD, level=c("SD", "Non-SD"))
df$LargeSD <- factor(df$LargeSD, level=c("LargeSD", "Non-LargeSD"))


###INV BRK###
categorical_boxplot(df$InvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nBreakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.brkpnt_10KSNP.png")
###LARGE INV BRK###
categorical_boxplot(df$LargeInvBrk, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge Inv Breakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Proportion", "happrop.vs.largebrkpnt_10KSNP.png")
###SD###
categorical_boxplot(df$SD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nSD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.SD_10KSNP.png")
###Large SD###
categorical_boxplot(df$LargeSD, df$Proportion_of_Unique_Haplotypes, "Proportion of Unique Haplotypes vs.\nLarge SD or Non-SD Window", "Window Type", "Haplotype Proportion", "happrop.vs.largeSD_10KSNP.png")



