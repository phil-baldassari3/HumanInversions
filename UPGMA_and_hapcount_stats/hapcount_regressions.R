library(ggplot2)
library(data.table)
library(car)
library(tidyverse)


get_p_from_lm <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#Function performs linear regression between 2 variables using Pearson's Correlation and a simple linear model y ~ x. Result is plotted and saved as png.
linear_regression_plot <- function(xdata, ydata, plottitle, xlabel, ylabel, outfile) {
  
  #Pearson's Correlation and linear model
  model <- lm(ydata ~ xdata)
  #prscor <- cor(ydata, xdata, method="pearson")
  
  anova_on_model <- anova(model)
  anova_df <- data.frame(anova_on_model)
  anova_df <- anova_df %>%
    select(Sum.Sq)
  summedsumsq <- sum(anova_df$Sum.Sq)
  anova_df$PEV <- anova_df$Sum.Sq / summedsumsq
  print(anova_df)
  
  #making ploting df
  df <- data.frame(xcol = xdata, ycol = ydata)
  
  #plotting
  res = 300
  scale = 1
  png(outfile, height = 3 * scale * res, width = 4 * scale * res, res = res)
  p <- ggplot(df, aes(x=xcol, y=ycol)) +
    geom_point(size=1, alpha=0.03) +
    geom_smooth(method='lm') + 
    #geom_text(x = 0.3, y = 0.4, size = 6, label = paste("Pearson's Correlation:", round(prscor, 3))) +
    #geom_text(x = 0.3, y = 0.36, size = 6, label = paste("R-squared:", round(summary(model)$adj.r.squared, 3))) +
    #geom_text(x = 0.3, y = 0.32, size = 6, label = paste("p-value:", round(get_p_from_lm(model), 3))) +
    labs(x=xlabel, y=ylabel, title=plottitle) +
    theme_classic()
  print(p)
  dev.off()
  
  
}


#Function performs linear regression between 2 variables using Pearson's Correlation and a simple linear model y ~ x. Result is plotted and saved as png.
categorical_regression_boxplot <- function(xdata, ydata, plottitle, xlabel, ylabel, outfile) {
  
  #linear model
  model <- lm(ydata ~ xdata)
  anova_on_model <- anova(model)
  anova_df <- data.frame(anova_on_model)
  anova_df <- anova_df %>%
    select(Sum.Sq)
  summedsumsq <- sum(anova_df$Sum.Sq)
  anova_df$PEV <- anova_df$Sum.Sq / summedsumsq
  print(anova_df)
  
  #making ploting df
  df <- data.frame(xcol = xdata, ycol = ydata)
  
  #plotting
  res = 300
  scale = 1.3
  png(outfile, height = 3 * scale * res, width = 4 * scale * res, res = res)
  p <- ggplot(df, aes(x=xcol, y=ycol)) +
    geom_boxplot(outliers = FALSE) +
    geom_smooth(method='lm', aes(group = 1)) + 
    labs(x=xlabel, y=ylabel, title=plottitle) +
    theme_classic()
  print(p)
  dev.off()
  
  res = 300
  scale = 1.3
  png(paste("density_", outfile), height = 3 * scale * res, width = 4 * scale * res, res = res)
  p2 <- ggplot(df, aes(x=ycol, fill=xcol)) +
    geom_density(position = "identity", alpha = 0.5, linewidth = 0.1) + #bins = 30
    scale_fill_discrete(name = "Window type") +
    labs(x=ylabel, title=plottitle) +
    theme_classic()
  print(p2)
  dev.off()
  
  
}



setwd("~/Desktop/Phadnis_Lab/human_inversions/uniq_hap_counts")

#"hapcount_SNPwindow/cleaned_autosomes_SNPwindow1000_SNPstep500_hapcount_snpden_labeled.tsv" "hapcount_BPwindow/cleaned_autosomes_BPwindow3000_BPstep1500_hapcount_snpden_labeled.tsv"
df <- as.data.frame(fread("hapcount_BPwindow/cleaned_autosomes_BPwindow3000_BPstep1500_hapcount_snpden_labeled.tsv",, sep="\t", header=TRUE))
df$inversion <- as.factor(df$inversion)
df$large_inversion <- as.factor(df$large_inversion)
df$inv_breakpoint <- as.factor(df$inv_breakpoint)
df$large_inv_breakpoint <- as.factor(df$large_inv_breakpoint)
head(df)





linear_regression_plot(df$NUM_uniq_haps, df$SNP_density, "SNP Density vs. Haplotype Count", "Haplotype Count", "SNPs/bp", "snpden.vs.branch.png")


categorical_regression_boxplot(df$inversion, df$NUM_uniq_haps, "Hap Count vs.\nInversion or Non-Inversion Window", "Window Type", "Haplotype Count", "branch.vs.inv.png")
categorical_regression_boxplot(df$large_inversion, df$NUM_uniq_haps, "Hap Count vs.\nInversion > 1Mb or Non-Inversion Window", "Window Type", "Haplotype Count", "branch.vs.largeinv.png")
categorical_regression_boxplot(df$inv_breakpoint, df$NUM_uniq_haps, "Hap Count vs.\nBreakpoint or Non-Breakpoint Window", "Window Type", "Haplotype Count", "branch.vs.brkpnt.png")
categorical_regression_boxplot(df$large_inv_breakpoint, df$NUM_uniq_haps, "Hap Count vs.\nBreakpoint (Inv > 1Mb) or Non-Breakpoint Window", "Window Type", "Haplotype Count", "branch.vs.largebrkpnt.png")









 # df$logitrepeat = logit(df$Repeat_Density)
# 
# pct75 <- quantile(df$NUM_uniq_haps, probs=0.75)
# pct90 <- quantile(df$NUM_uniq_haps, probs=0.9)
# pct98 <- quantile(df$NUM_uniq_haps, probs=0.98)
# 
# df75 <- df[df$NUM_uniq_haps > pct75,]
# df90 <- df[df$NUM_uniq_haps > pct90,]
# df98 <- df[df$NUM_uniq_haps > pct98,]
# 
# linear_regression_plot(df$NUM_uniq_haps, df$SNP_Density, "SNP Density vs. Branch Length", "Haplotype Count", "SNPs/bp", "snpden_vs_avgbranchlen.png")
# linear_regression_plot(df75$NUM_uniq_haps, df75$SNP_Density, "SNP Density vs. Branch Length for windows in 75th percentile of Branch Length", "Haplotype Count", "SNPs/bp", "snpden_vs_avgbranchlen_75.png")
# linear_regression_plot(df90$NUM_uniq_haps, df90$SNP_Density, "SNP Density vs. Branch Length for windows in 90th percentile of Branch Length", "Haplotype Count", "SNPs/bp", "snpden_vs_avgbranchlen_90.png")
# linear_regression_plot(df98$NUM_uniq_haps, df98$SNP_Density, "SNP Density vs. Branch Length for windows in 98th percentile of Branch Length", "Haplotype Count", "SNPs/bp", "snpden_vs_avgbranchlen_98.png")
# 
# linear_regression_plot(df$NUM_uniq_haps, df$Repeat_Density, "Repeat Density vs. Branch Length", "Haplotype Count", "Repeat Density", "repeatden_vs_avgbranchlen.png")
# linear_regression_plot(df75$NUM_uniq_haps, df75$Repeat_Density, "Repeat Density vs. Branch Length for windows in 75th percentile of Branch Length", "Haplotype Count", "Repeat Density", "repeatden_vs_avgbranchlen_75.png")
# linear_regression_plot(df90$NUM_uniq_haps, df90$Repeat_Density, "Repeat Density vs. Branch Length for windows in 90th percentile of Branch Length", "Haplotype Count", "Repeat Density", "repeatden_vs_avgbranchlen_90.png")
# linear_regression_plot(df98$NUM_uniq_haps, df98$Repeat_Density, "Repeat Density vs. Branch Length for windows in 98th percentile of Branch Length", "Haplotype Count", "Repeat Density", "repeatden_vs_avgbranchlen_98.png")

#linear_regression_plot(df$NUM_uniq_haps, df$logitrepeat, "Repeat Density vs. Branch Length", "Haplotype Count", "Repeat Density (Logit Transformed)", "repeatden_vs_avgbranchlen_logit.png")

