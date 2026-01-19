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



setwd("/Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/hierarchical_clustering/some_unfinished_100SNPwin_csvs")

df <- as.data.frame(fread("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden_labeled.tsv", sep="\t", header=TRUE))
df$inversion <- as.factor(df$inversion)
df$large_inversion <- as.factor(df$large_inversion)
df$inv_breakpoint <- as.factor(df$inv_breakpoint)
df$large_inv_breakpoint <- as.factor(df$large_inv_breakpoint)
head(df)





linear_regression_plot(df$Avg_Branch_Length, df$SNP_Density, "SNP Density vs. Branch Length", "Average Branch Length", "SNPs/bp", "snpden.vs.branch.png")
linear_regression_plot(df$Avg_Branch_Length, df$Repeat_Density, "Repeat Density vs. Branch Length", "Average Branch Length", "Repeat Density", "repeatden.vs.branch.png")
linear_regression_plot(df$SNP_Density, df$Repeat_Density, "Repeat Density vs. SNP Density", "SNPs/bp", "Repeat Density", "repeatden.vs.snpden.png")


categorical_regression_boxplot(df$inversion, df$Avg_Branch_Length, "Avg Branch Length vs.\nInversion or Non-Inversion Window", "Window Type", "Average Branch Length", "branch.vs.inv.png")
categorical_regression_boxplot(df$large_inversion, df$Avg_Branch_Length, "Avg Branch Length vs.\nInversion > 1Mb or Non-Inversion Window", "Window Type", "Average Branch Length", "branch.vs.largeinv.png")
categorical_regression_boxplot(df$inv_breakpoint, df$Avg_Branch_Length, "Avg Branch Length vs.\nBreakpoint or Non-Breakpoint Window", "Window Type", "Average Branch Length", "branch.vs.brkpnt.png")
categorical_regression_boxplot(df$large_inv_breakpoint, df$Avg_Branch_Length, "Avg Branch Length vs.\nBreakpoint (Inv > 1Mb) or Non-Breakpoint Window", "Window Type", "Average Branch Length", "branch.vs.largebrkpnt.png")









 # df$logitrepeat = logit(df$Repeat_Density)
# 
# pct75 <- quantile(df$Avg_Branch_Length, probs=0.75)
# pct90 <- quantile(df$Avg_Branch_Length, probs=0.9)
# pct98 <- quantile(df$Avg_Branch_Length, probs=0.98)
# 
# df75 <- df[df$Avg_Branch_Length > pct75,]
# df90 <- df[df$Avg_Branch_Length > pct90,]
# df98 <- df[df$Avg_Branch_Length > pct98,]
# 
# linear_regression_plot(df$Avg_Branch_Length, df$SNP_Density, "SNP Density vs. Branch Length", "Average Branch Length", "SNPs/bp", "snpden_vs_avgbranchlen.png")
# linear_regression_plot(df75$Avg_Branch_Length, df75$SNP_Density, "SNP Density vs. Branch Length for windows in 75th percentile of Branch Length", "Average Branch Length", "SNPs/bp", "snpden_vs_avgbranchlen_75.png")
# linear_regression_plot(df90$Avg_Branch_Length, df90$SNP_Density, "SNP Density vs. Branch Length for windows in 90th percentile of Branch Length", "Average Branch Length", "SNPs/bp", "snpden_vs_avgbranchlen_90.png")
# linear_regression_plot(df98$Avg_Branch_Length, df98$SNP_Density, "SNP Density vs. Branch Length for windows in 98th percentile of Branch Length", "Average Branch Length", "SNPs/bp", "snpden_vs_avgbranchlen_98.png")
# 
# linear_regression_plot(df$Avg_Branch_Length, df$Repeat_Density, "Repeat Density vs. Branch Length", "Average Branch Length", "Repeat Density", "repeatden_vs_avgbranchlen.png")
# linear_regression_plot(df75$Avg_Branch_Length, df75$Repeat_Density, "Repeat Density vs. Branch Length for windows in 75th percentile of Branch Length", "Average Branch Length", "Repeat Density", "repeatden_vs_avgbranchlen_75.png")
# linear_regression_plot(df90$Avg_Branch_Length, df90$Repeat_Density, "Repeat Density vs. Branch Length for windows in 90th percentile of Branch Length", "Average Branch Length", "Repeat Density", "repeatden_vs_avgbranchlen_90.png")
# linear_regression_plot(df98$Avg_Branch_Length, df98$Repeat_Density, "Repeat Density vs. Branch Length for windows in 98th percentile of Branch Length", "Average Branch Length", "Repeat Density", "repeatden_vs_avgbranchlen_98.png")

#linear_regression_plot(df$Avg_Branch_Length, df$logitrepeat, "Repeat Density vs. Branch Length", "Average Branch Length", "Repeat Density (Logit Transformed)", "repeatden_vs_avgbranchlen_logit.png")

