library(ggplot2)
library(forcats)


setwd("~/Desktop/Phadnis_Lab/human_inversions/testing_hapcount_classification_power/bootstrap_LD")

means_dprime_df <- read.csv("means_dprime_comparisons_with_all_samples.csv")


### GEOMETRIC MEAN ####

plot_dprime_gmeans_bargraph <- function(df, outpng){
    ###Function plots a bar graph of Gmean(\D'\) for each comparison 
    ###ordered from greatest to least and colored by category
    
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1.5
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=fct_reorder(comparison_name, -dprime_gmean), y=dprime_gmean, fill=comparison_type)) +
        geom_col() +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Geometric Mean of |D'| by Comparison") + xlab("Comparison") + ylab("Geometric_Mean(|D'|)") +
        theme_bw() + theme(axis.text.x = element_blank()) #, axis.title.x = element_blank())
    print(p)
    dev.off()
    
}


plot_gmean_dprime_by_category <- function(df, outpng){
    ###Function plots boxplots of Gmean(\D'\) distributions for each category
    
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    model <- lm(data=df, dprime_gmean ~ comparison_type)
    print(summary(model))
    print(anova(model))
    
    res = 300
    scale = 1.5
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=comparison_type, y=dprime_gmean, fill=comparison_type)) +
        geom_boxplot() +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Geometric Mean of |D'| by Comparison Type") + xlab("Comparison Type") + ylab("Geometric_Mean(|D'|)") +
        theme_bw() + theme(legend.position = "none")
    print(p)
    dev.off()
}


plot_gmeanDprime_to_distance_regression <- function(df, outpng){
    ###Function plots a linear regression of comparison Gmean(\D'\) vs. distance between peaks
    
    #df <- df[df$dprime_gmean < 0.25, ]
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    df$a_midpoint <- (df$endA - df$startA) / 2
    df$b_midpoint <- (df$endB - df$startB) / 2
    df$Distance <- abs(floor((df$b_midpoint - df$a_midpoint) / 2))
    
    model <- lm(data=df, dprime_gmean ~ Distance)
    print(summary(model))
    
    res = 300
    scale = 1.5
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=Distance, y=dprime_gmean, color=comparison_type)) +
        geom_point(size=2, alpha=0.7) + geom_smooth(method="lm", se=FALSE, linewidth=0.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Regression of Gmean of |D'| to Peak Distance") + xlab("Distance b/n Peaks") + ylab("Geometric_Mean(|D'|)") +
        theme_bw()
    print(p)
    dev.off()
    
    
}


plot_gmeanDprime_to_sitecomps_regression <- function(df, outpng){
    ###Function plots a linear regression of comparison Gmean(\D'\) vs. number of site comparisons
    
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    model <- lm(data=df, dprime_gmean ~ Num_Site_Comparisons)
    print(summary(model))
    
    res = 300
    scale = 1.5
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=Num_Site_Comparisons, y=dprime_gmean, color=comparison_type)) +
        geom_point(size=2, alpha=0.7) + geom_smooth(method="lm", se=FALSE, linewidth=0.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Regression of Gmean of |D'| to # of Site Comparisons") + xlab("# of Site Comparisons") + ylab("Geometric_Mean(|D'|)") +
        theme_bw()
    print(p)
    dev.off()
}








plot_dprime_gmeans_bargraph(means_dprime_df, "GmeanDprime_by_comparison.png")

plot_gmean_dprime_by_category(means_dprime_df, "GmeanDprime_by_comparison_type.png")

plot_gmeanDprime_to_distance_regression(means_dprime_df, "Regression_GmeanDprime_vs_peakDistance.png")
  
plot_gmeanDprime_to_sitecomps_regression(means_dprime_df, "Regression_GmeanDprime_vs_num_site_comps.png")




















###########################################################################################################################

### MEAN ###


plot_dprime_means_bargraph <- function(df, outpng){
    ###Function plots a bar graph of Gmean(\D'\) for each comparison 
    ###ordered from greatest to least and colored by category
    
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1.5
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=fct_reorder(comparison_name, -dprime_mean), y=dprime_mean, fill=comparison_type)) +
        geom_col() +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Mean of |D'| by Comparison") + xlab("Comparison") + ylab("Mean(|D'|)") +
        theme_bw() + theme(axis.text.x = element_blank()) #, axis.title.x = element_blank())
    print(p)
    dev.off()
    
}


plot_mean_dprime_by_category <- function(df, outpng){
    ###Function plots boxplots of Gmean(\D'\) distributions for each category
    
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    model <- lm(data=df, dprime_mean ~ comparison_type)
    print(summary(model))
    print(anova(model))
    
    res = 300
    scale = 1.5
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=comparison_type, y=dprime_mean, fill=comparison_type)) +
        geom_boxplot() +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Mean of |D'| by Comparison Type") + xlab("Comparison Type") + ylab("Mean(|D'|)") +
        theme_bw() + theme(legend.position = "none")
    print(p)
    dev.off()
}


plot_meanDprime_to_distance_regression <- function(df, outpng){
    ###Function plots a linear regression of comparison Gmean(\D'\) vs. distance between peaks
    
    #df <- df[df$dprime_mean < 0.25, ]
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    df$a_midpoint <- (df$endA - df$startA) / 2
    df$b_midpoint <- (df$endB - df$startB) / 2
    df$Distance <- abs(floor((df$b_midpoint - df$a_midpoint) / 2))
    
    model <- lm(data=df, dprime_mean ~ Distance)
    print(summary(model))
    
    res = 300
    scale = 1.5
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=Distance, y=dprime_mean, color=comparison_type)) +
        geom_point(size=2, alpha=0.7) + geom_smooth(method="lm", se=FALSE, linewidth=0.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Regression of Mean of |D'| to Peak Distance") + xlab("Distance b/n Peaks") + ylab("Mean(|D'|)") +
        theme_bw()
    print(p)
    dev.off()
    
    
}


plot_meanDprime_to_sitecomps_regression <- function(df, outpng){
    ###Function plots a linear regression of comparison Gmean(\D'\) vs. number of site comparisons
    
    df <- df[df$comparison_type != "Adjacent INVs", ]
    
    model <- lm(data=df, dprime_mean ~ Num_Site_Comparisons)
    print(summary(model))
    
    res = 300
    scale = 1.5
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=Num_Site_Comparisons, y=dprime_mean, color=comparison_type)) +
        geom_point(size=2, alpha=0.7) + geom_smooth(method="lm", se=FALSE, linewidth=0.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        ggtitle("Regression of Mean of |D'| to # of Site Comparisons") + xlab("# of Site Comparisons") + ylab("Mean(|D'|)") +
        theme_bw()
    print(p)
    dev.off()
}






plot_dprime_means_bargraph(means_dprime_df, "MeanDprime_by_comparison.png")

plot_mean_dprime_by_category(means_dprime_df, "MeanDprime_by_comparison_type.png")

plot_meanDprime_to_distance_regression(means_dprime_df, "Regression_MeanDprime_vs_peakDistance.png")

plot_meanDprime_to_sitecomps_regression(means_dprime_df, "Regression_MeanDprime_vs_num_site_comps.png")

