library(ggplot2)
library(ggridges)
library(forcats)
library(RColorBrewer)
library(tidytext)


setwd("~/Desktop/Phadnis_Lab/human_inversions/testing_hapcount_classification_power/bootstrap_LD")

df <- read.csv("bootstrapLD_all_comparisions.csv")
df2 <- read.csv("individual_bootstapped_LDstats.csv")



### GEOMETRIC MEAN ####

plot_Dprime_by_type <- function(bootstrap_df, outpng){
    ### Function plots a faceted plot of Dprime distributions by comparision type and saves an output png
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=bootstrap_df, aes(x=dprime_gmean, fill=comparison_type)) +
        geom_histogram(bins=30, aes(y = after_stat(density))) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Bootstrap Gmean(D')") + ylab("Density") + 
        ggtitle("Distribution of Bootstrap D' by Comparison Type") +
        theme_bw() +
        facet_wrap(~comparison_type, ncol=1)
    print(p)
    dev.off()
    
}


ttest_Dprime_by_type <- function(bootstrap_df){
    ### Runs a Welch's T-test on Dprimes between INV to INV vs. other types
    
    #filtering dfs
    df_inv <- bootstrap_df[bootstrap_df$comparison_type == "INV to INV", ]
    df_other <- bootstrap_df[bootstrap_df$comparison_type != "INV to INV", ]
    
    #grabbing vectors
    inv_dprimes <- df_inv$dprime_gmean
    other_dprimes <- df_other$dprime_gmean
    
    #t-test
    t <- t.test(inv_dprimes, other_dprimes, var.equal=FALSE)
    print("Is there a significant difference in means between Dprime of INV to INV vs Other types")
    print(t)

    
}



plot_Dprime_by_comparison <- function(bootstrap_df, outpng){
    ### Function plots a Dprime boxplots and distributions by comparision name and saves an output png
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(paste0("box_", outpng), height = 3 * scale * res, width = 4 * scale * res, res = res)
    p1 <- ggplot(data=bootstrap_df, aes(x=fct_reorder(comparison_name, -dprime_gmean, .fun=median), y=dprime_gmean, fill=comparison_type)) +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Comparision") + ylab("Bootstrap Gmean(D')") + 
        ggtitle("Bootstrap D' by Comparison") +
        theme_bw() + theme(axis.text.x = element_blank())
    print(p1)
    dev.off()
    
    res = 300
    scale = 3
    png(paste0("dist_", outpng), height = 4 * scale * res, width = 3 * scale * res, res = res)
    p2 <- ggplot(data=bootstrap_df, aes(x=dprime_gmean, y=fct_reorder(comparison_name, -dprime_gmean, .fun=median), height=after_stat(density), fill=comparison_type)) +
        geom_density_ridges(stat = "binline", bins = 30, scale = 1.5, draw_baseline = FALSE) +
        #geom_density_ridges(stat = "density", scale = 1.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        ylab("Comparision") + xlab("Bootstrap Gmean(D')") + 
        ggtitle("Bootstrap D' by Comparison") +
        theme_bw()
    print(p2)
    dev.off()
    
}


plot_IndvLD_by_type <- function(indv_df, outpng){
    ### Function plots a faceted plot of Individual LDstat_gmean distributions by comparision type and saves an output png
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=LDstat_gmean, fill=comparison_type)) +
        geom_histogram(bins=30, aes(y = after_stat(density))) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Individual LDstat (from gmean)") + ylab("Density") + 
        ggtitle("Distribution of Individual LDstats by Comparison Type") +
        theme_bw() +
        facet_wrap(~comparison_type, ncol=1)
    print(p)
    dev.off()
    
}



ttest_IndvLD_by_type <- function(indv_df){
    ### Runs a Welch's T-test on LDstat_gmeans between INV to INV vs. other types
    
    #filtering dfs
    df_inv <- indv_df[indv_df$comparison_type == "INV to INV", ]
    df_other <- indv_df[indv_df$comparison_type != "INV to INV", ]
    
    #grabbing vectors
    inv_LD <- df_inv$LDstat_gmean
    other_LD <- df_other$LDstat_gmean
    
    #t-test
    t <- t.test(inv_LD, other_LD, var.equal=FALSE)
    print("Is there a significant difference in means between LDstats of INV to INV vs Other types")
    print(t)
    
    
}



plot_IndvLD_by_comparison <- function(indv_df, outpng){
    ### Function plots a Dprime boxplots and distributions by comparision name and saves an output png
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(paste0("box_", outpng), height = 3 * scale * res, width = 4 * scale * res, res = res)
    p1 <- ggplot(data=indv_df, aes(x=fct_reorder(comparison_name, -LDstat_gmean, .fun=median), y=LDstat_gmean, fill=comparison_type)) +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Comparision") + ylab("Individual LDstat (from gmean)") + 
        ggtitle("LDstats by Comparison") +
        theme_bw() + theme(axis.text.x = element_blank())
    print(p1)
    dev.off()
    
    res = 300
    scale = 3
    png(paste0("dist_", outpng), height = 4 * scale * res, width = 3 * scale * res, res = res)
    p2 <- ggplot(data=indv_df, aes(x=LDstat_gmean, y=fct_reorder(comparison_name, -LDstat_gmean, .fun=median), height=stat(density), fill=comparison_type)) +
        geom_density_ridges(stat = "binline", bins = 30, scale = 1.5, draw_baseline = FALSE) +
        #geom_density_ridges(stat = "density", scale = 1.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        ylab("Comparision") + xlab("Individual LDstat (from gmean)") + 
        ggtitle("LDstats by Comparison") +
        theme_bw()
    print(p2)
    dev.off()
    
}



plot_IndvLD_by_comparison_experimental <- function(indv_df, outpng){
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1
    png(outpng, height = 9 * scale * res, width = 16 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=reorder_within(Sample, LDstat_gmean, comparison_name), y=LDstat_gmean, color=comparison_type)) +
        geom_point(size=1.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(color = "Comparison Type") +
        xlab("Individual") + ylab("Individual LDstat (from gmean)") + 
        ggtitle("LDstats by Comparison") +
        facet_wrap(~comparison_name, ncol=9, scales="free_x") +
        theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank())
    print(p)
    dev.off()
    
}

plot_IndvLD_ecdf_by_comparison <- function(indv_df, outpng){
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1
    png(outpng, height = 9 * scale * res, width = 16 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=LDstat_gmean, color=comparison_type)) +
        stat_ecdf(geom = "step", linewidth=1.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(color = "Comparison Type") +
        xlab("Individual LDstat (from gmean)") + ylab("ecdf") + 
        ggtitle("eCDF of LDstats by Comparison") +
        facet_wrap(~comparison_name, ncol=9) +
        theme_bw()
    print(p)
    dev.off()
    
}




regression_dprimes_to_distance <- function(bootstrap_df, outpng){
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    #compute distances
    bootstrap_df$midA <- bootstrap_df$endA - bootstrap_df$startA
    bootstrap_df$midB <- bootstrap_df$endB - bootstrap_df$startB
    bootstrap_df$distance <- abs(bootstrap_df$midB - bootstrap_df$midA)
    
    #linear model
    model <- lm(dprime_gmean ~ distance, data=bootstrap_df)
    print(summary(model))
    
    #plot
    p <- ggplot(data=bootstrap_df, aes(x=distance, y=dprime_gmean, color=comparison_type)) +
        geom_point(alpha=0.5) +
        geom_smooth(method = "lm") +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        xlab("Distance between peaks (bp)") + ylab("Bootstrap D'") + ggtitle("Regression of Bootstrap D' to Distance b/n Peaks") +
        theme_bw()
    print(p)
    
}


regression_dprimes_to_numsitecomps <- function(bootstrap_df, outpng){
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    #linear model
    model <- lm(dprime_gmean ~ Num_Site_Comparisons, data=bootstrap_df)
    print(summary(model))
    anva <- anova(model)
    print(anva)
    
    #plot
    p <- ggplot(data=bootstrap_df, aes(x=Num_Site_Comparisons, y=dprime_gmean, color=comparison_type)) +
        geom_point(alpha=0.5) +
        geom_smooth(method = "lm") +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        xlab("Number of Site Comparisions") + ylab("Bootstrap D'") + ggtitle("Regression of Bootstrap D' to Number of Site Comparisons") +
        theme_bw()
    print(p)
    
}


regression_dprimes_to_numsitecomps_facet <- function(bootstrap_df, outpng){
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    #linear model
    model <- lm(dprime_gmean ~ Num_Site_Comparisons, data=bootstrap_df)
    print(summary(model))
    
    #plot
    p <- ggplot(data=bootstrap_df, aes(x=Num_Site_Comparisons, y=dprime_gmean, color=comparison_type)) +
        geom_point(alpha=0.5) +
        geom_smooth(method = "lm") +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        xlab("Number of Site Comparisions") + ylab("Bootstrap D'") + ggtitle("Regression of Bootstrap D' to Number of Site Comparisons") +
        facet_wrap(~comparison_name, ncol=9) +
        theme_bw()
        
    print(p)
    
}






plot_Dprime_by_type(df, "bootstrap_Dprimes_gmean_by_type.png")

ttest_Dprime_by_type(df)

plot_Dprime_by_comparison(df, "bootstrap_Dprimes_gmean_by_comparison.png")

plot_IndvLD_by_type(df2, "individual_LDstats_gmean_by_type.png")

ttest_IndvLD_by_type(df2)

plot_IndvLD_by_comparison(df2, "individual_LDstats_gmean_by_type.png")

plot_IndvLD_by_comparison_experimental(df2, "individual_LDstats_gmean_sorted.png")

plot_IndvLD_ecdf_by_comparison(df2, "individual_LDstats_gmean_ecdf.png")




# regression_dprimes_to_distance(df, "")
# regression_dprimes_to_numsitecomps(df, "")
# regression_dprimes_to_numsitecomps_facet(df, "")




###########################################################################################################################


### MEAN ####

plot_Dprime_by_type <- function(bootstrap_df, outpng){
    ### Function plots a faceted plot of Dprime distributions by comparision type and saves an output png
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=bootstrap_df, aes(x=dprime_mean, fill=comparison_type)) +
        geom_histogram(bins=30, aes(y = after_stat(density))) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Bootstrap Mean(D')") + ylab("Density") + 
        ggtitle("Distribution of Bootstrap D' by Comparison Type") +
        theme_bw() +
        facet_wrap(~comparison_type, ncol=1)
    print(p)
    dev.off()
    
}


ttest_Dprime_by_type <- function(bootstrap_df){
    ### Runs a Welch's T-test on Dprimes between INV to INV vs. other types
    
    #filtering dfs
    df_inv <- bootstrap_df[bootstrap_df$comparison_type == "INV to INV", ]
    df_other <- bootstrap_df[bootstrap_df$comparison_type != "INV to INV", ]
    
    #grabbing vectors
    inv_dprimes <- df_inv$dprime_mean
    other_dprimes <- df_other$dprime_mean
    
    #t-test
    t <- t.test(inv_dprimes, other_dprimes, var.equal=FALSE)
    print("Is there a significant difference in means between Dprime of INV to INV vs Other types")
    print(t)
    
    
}



plot_Dprime_by_comparison <- function(bootstrap_df, outpng){
    ### Function plots a Dprime boxplots and distributions by comparision name and saves an output png
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(paste0("box_", outpng), height = 3 * scale * res, width = 4 * scale * res, res = res)
    p1 <- ggplot(data=bootstrap_df, aes(x=fct_reorder(comparison_name, -dprime_mean, .fun=median), y=dprime_mean, fill=comparison_type)) +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Comparision") + ylab("Bootstrap Mean(D')") + 
        ggtitle("Bootstrap D' by Comparison") +
        theme_bw() + theme(axis.text.x = element_blank())
    print(p1)
    dev.off()
    
    res = 300
    scale = 3
    png(paste0("dist_", outpng), height = 4 * scale * res, width = 3 * scale * res, res = res)
    p2 <- ggplot(data=bootstrap_df, aes(x=dprime_mean, y=fct_reorder(comparison_name, -dprime_mean, .fun=median), height=after_stat(density), fill=comparison_type)) +
        geom_density_ridges(stat = "binline", bins = 30, scale = 1.5, draw_baseline = FALSE) +
        #geom_density_ridges(stat = "density", scale = 1.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        ylab("Comparision") + xlab("Bootstrap Mean(D')") + 
        ggtitle("Bootstrap D' by Comparison") +
        theme_bw()
    print(p2)
    dev.off()
    
}


plot_IndvLD_by_type <- function(indv_df, outpng){
    ### Function plots a faceted plot of Individual LDstat_mean distributions by comparision type and saves an output png
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=LDstat_mean, fill=comparison_type)) +
        geom_histogram(bins=30, aes(y = after_stat(density))) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Individual LDstat (from mean)") + ylab("Density") + 
        ggtitle("Distribution of Individual LDstats by Comparison Type") +
        theme_bw() +
        facet_wrap(~comparison_type, ncol=1)
    print(p)
    dev.off()
    
}



ttest_IndvLD_by_type <- function(indv_df){
    ### Runs a Welch's T-test on LDstat_means between INV to INV vs. other types
    
    #filtering dfs
    df_inv <- indv_df[indv_df$comparison_type == "INV to INV", ]
    df_other <- indv_df[indv_df$comparison_type != "INV to INV", ]
    
    #grabbing vectors
    inv_LD <- df_inv$LDstat_mean
    other_LD <- df_other$LDstat_mean
    
    #t-test
    t <- t.test(inv_LD, other_LD, var.equal=FALSE)
    print("Is there a significant difference in means between LDstat_means of INV to INV vs Other types")
    print(t)
    
    
}



plot_IndvLD_by_comparison <- function(indv_df, outpng){
    ### Function plots a Dprime boxplots and distributions by comparision name and saves an output png
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(paste0("box_", outpng), height = 3 * scale * res, width = 4 * scale * res, res = res)
    p1 <- ggplot(data=indv_df, aes(x=fct_reorder(comparison_name, -LDstat_mean, .fun=median), y=LDstat_mean, fill=comparison_type)) +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Comparision") + ylab("Individual LDstat (from mean)") + 
        ggtitle("LDstats by Comparison") +
        theme_bw() + theme(axis.text.x = element_blank())
    print(p1)
    dev.off()
    
    res = 300
    scale = 3
    png(paste0("dist_", outpng), height = 4 * scale * res, width = 3 * scale * res, res = res)
    p2 <- ggplot(data=indv_df, aes(x=LDstat_mean, y=fct_reorder(comparison_name, -LDstat_mean, .fun=median), height=stat(density), fill=comparison_type)) +
        geom_density_ridges(stat = "binline", bins = 30, scale = 1.5, draw_baseline = FALSE) +
        #geom_density_ridges(stat = "density", scale = 1.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        ylab("Comparision") + xlab("Individual LDstat (from mean)") + 
        ggtitle("LDstats by Comparison") +
        theme_bw()
    print(p2)
    dev.off()
    
}



plot_IndvLD_by_comparison_experimental <- function(indv_df, outpng){
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1
    png(outpng, height = 9 * scale * res, width = 16 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=reorder_within(Sample, LDstat_mean, comparison_name), y=LDstat_mean, color=comparison_type)) +
        geom_point(size=1.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(color = "Comparison Type") +
        xlab("Individual") + ylab("Individual LDstat (from mean)") + 
        ggtitle("LDstats by Comparison") +
        facet_wrap(~comparison_name, ncol=9, scales="free_x") +
        theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank())
    print(p)
    dev.off()
    
}

plot_IndvLD_ecdf_by_comparison <- function(indv_df, outpng){
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1
    png(outpng, height = 9 * scale * res, width = 16 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=LDstat_mean, color=comparison_type)) +
        stat_ecdf(geom = "step", linewidth=1.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(color = "Comparison Type") +
        xlab("Individual LDstat (from mean)") + ylab("ecdf") + 
        ggtitle("eCDF of LDstats by Comparison") +
        facet_wrap(~comparison_name, ncol=9) +
        theme_bw()
    print(p)
    dev.off()
    
}






plot_Dprime_by_type(df, "bootstrap_Dprimes_mean_by_type.png")

ttest_Dprime_by_type(df)

plot_Dprime_by_comparison(df, "bootstrap_Dprimes_mean_by_comparison.png")

plot_IndvLD_by_type(df2, "individual_LDstats_mean_by_type.png")

ttest_IndvLD_by_type(df2)

plot_IndvLD_by_comparison(df2, "individual_LDstats_mean_by_type.png")

plot_IndvLD_by_comparison_experimental(df2, "individual_LDstats_mean_sorted.png")

plot_IndvLD_ecdf_by_comparison(df2, "individual_LDstats_mean_ecdf.png")







###########################################################################################################################


### MEAN NON-ZERO ####

plot_Dprime_by_type <- function(bootstrap_df, outpng){
    ### Function plots a faceted plot of Dprime distributions by comparision type and saves an output png
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=bootstrap_df, aes(x=dprime_mean_nonzero, fill=comparison_type)) +
        geom_histogram(bins=30, aes(y = after_stat(density))) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Bootstrap Mean(D'>0)") + ylab("Density") + 
        ggtitle("Distribution of Bootstrap D' by Comparison Type") +
        theme_bw() +
        facet_wrap(~comparison_type, ncol=1)
    print(p)
    dev.off()
    
}


ttest_Dprime_by_type <- function(bootstrap_df){
    ### Runs a Welch's T-test on Dprimes between INV to INV vs. other types
    
    #filtering dfs
    df_inv <- bootstrap_df[bootstrap_df$comparison_type == "INV to INV", ]
    df_other <- bootstrap_df[bootstrap_df$comparison_type != "INV to INV", ]
    
    #grabbing vectors
    inv_dprimes <- df_inv$dprime_mean_nonzero
    other_dprimes <- df_other$dprime_mean_nonzero
    
    #t-test
    t <- t.test(inv_dprimes, other_dprimes, var.equal=FALSE)
    print("Is there a significant difference in means between Dprime of INV to INV vs Other types")
    print(t)
    
    
}



plot_Dprime_by_comparison <- function(bootstrap_df, outpng){
    ### Function plots a Dprime boxplots and distributions by comparision name and saves an output png
    
    bootstrap_df <- bootstrap_df[bootstrap_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(paste0("box_", outpng), height = 3 * scale * res, width = 4 * scale * res, res = res)
    p1 <- ggplot(data=bootstrap_df, aes(x=fct_reorder(comparison_name, -dprime_mean_nonzero, .fun=median), y=dprime_mean_nonzero, fill=comparison_type)) +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Comparision") + ylab("Bootstrap Mean(D'>0)") + 
        ggtitle("Bootstrap D' by Comparison") +
        theme_bw() + theme(axis.text.x = element_blank())
    print(p1)
    dev.off()
    
    res = 300
    scale = 3
    png(paste0("dist_", outpng), height = 4 * scale * res, width = 3 * scale * res, res = res)
    p2 <- ggplot(data=bootstrap_df, aes(x=dprime_mean_nonzero, y=fct_reorder(comparison_name, -dprime_mean_nonzero, .fun=median), height=after_stat(density), fill=comparison_type)) +
        geom_density_ridges(stat = "binline", bins = 30, scale = 1.5, draw_baseline = FALSE) +
        #geom_density_ridges(stat = "density", scale = 1.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        ylab("Comparision") + xlab("Bootstrap Mean(D'>0)") + 
        ggtitle("Bootstrap D' by Comparison") +
        theme_bw()
    print(p2)
    dev.off()
    
}


plot_IndvLD_by_type <- function(indv_df, outpng){
    ### Function plots a faceted plot of Individual LDstat_mean_nonzero distributions by comparision type and saves an output png
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(outpng, height = 4 * scale * res, width = 3 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=LDstat_mean_nonzero, fill=comparison_type)) +
        geom_histogram(bins=30, aes(y = after_stat(density))) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Individual LDstat (from mean nonzero)") + ylab("Density") + 
        ggtitle("Distribution of Individual LDstats by Comparison Type") +
        theme_bw() +
        facet_wrap(~comparison_type, ncol=1)
    print(p)
    dev.off()
    
}



ttest_IndvLD_by_type <- function(indv_df){
    ### Runs a Welch's T-test on LDstat_mean_nonzeros between INV to INV vs. other types
    
    #filtering dfs
    df_inv <- indv_df[indv_df$comparison_type == "INV to INV", ]
    df_other <- indv_df[indv_df$comparison_type != "INV to INV", ]
    
    #grabbing vectors
    inv_LD <- df_inv$LDstat_mean_nonzero
    other_LD <- df_other$LDstat_mean_nonzero
    
    #t-test
    t <- t.test(inv_LD, other_LD, var.equal=FALSE)
    print("Is there a significant difference in means between LDstat_mean_nonzeros of INV to INV vs Other types")
    print(t)
    
    
}



plot_IndvLD_by_comparison <- function(indv_df, outpng){
    ### Function plots a Dprime boxplots and distributions by comparision name and saves an output png
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 2
    png(paste0("box_", outpng), height = 3 * scale * res, width = 4 * scale * res, res = res)
    p1 <- ggplot(data=indv_df, aes(x=fct_reorder(comparison_name, -LDstat_mean_nonzero, .fun=median), y=LDstat_mean_nonzero, fill=comparison_type)) +
        geom_boxplot(outlier.size=0.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        xlab("Comparision") + ylab("Individual LDstat (from mean nonzero)") + 
        ggtitle("LDstats by Comparison") +
        theme_bw() + theme(axis.text.x = element_blank())
    print(p1)
    dev.off()
    
    res = 300
    scale = 3
    png(paste0("dist_", outpng), height = 4 * scale * res, width = 3 * scale * res, res = res)
    p2 <- ggplot(data=indv_df, aes(x=LDstat_mean_nonzero, y=fct_reorder(comparison_name, -LDstat_mean_nonzero, .fun=median), height=stat(density), fill=comparison_type)) +
        geom_density_ridges(stat = "binline", bins = 30, scale = 1.5, draw_baseline = FALSE) +
        #geom_density_ridges(stat = "density", scale = 1.5) +
        scale_fill_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(fill = "Comparison Type") +
        ylab("Comparision") + xlab("Individual LDstat (from mean nonzero)") + 
        ggtitle("LDstats by Comparison") +
        theme_bw()
    print(p2)
    dev.off()
    
}



plot_IndvLD_by_comparison_experimental <- function(indv_df, outpng){
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1
    png(outpng, height = 9 * scale * res, width = 16 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=reorder_within(Sample, LDstat_mean_nonzero, comparison_name), y=LDstat_mean_nonzero, color=comparison_type)) +
        geom_point(size=1.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(color = "Comparison Type") +
        xlab("Individual") + ylab("Individual LDstat (from mean nonzero)") + 
        ggtitle("LDstats by Comparison") +
        facet_wrap(~comparison_name, ncol=9, scales="free_x") +
        theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                           panel.grid.major.x = element_blank(),
                           panel.grid.minor.x = element_blank())
    print(p)
    dev.off()
    
}

plot_IndvLD_ecdf_by_comparison <- function(indv_df, outpng){
    
    indv_df <- indv_df[indv_df$comparison_type != "Adjacent INVs", ]
    
    res = 300
    scale = 1
    png(outpng, height = 9 * scale * res, width = 16 * scale * res, res = res)
    p <- ggplot(data=indv_df, aes(x=LDstat_mean_nonzero, color=comparison_type)) +
        stat_ecdf(geom = "step", linewidth=1.5) +
        scale_color_manual(values=c("#FC8D62", "#66C2A5", "#8DA0CB")) +
        labs(color = "Comparison Type") +
        xlab("Individual LDstat (from mean nonzero)") + ylab("ecdf") + 
        ggtitle("eCDF of LDstats by Comparison") +
        facet_wrap(~comparison_name, ncol=9) +
        theme_bw()
    print(p)
    dev.off()
    
}






plot_Dprime_by_type(df, "bootstrap_Dprimes_mean_nonzero_by_type.png")

ttest_Dprime_by_type(df)

plot_Dprime_by_comparison(df, "bootstrap_Dprimes_mean_nonzero_by_comparison.png")

plot_IndvLD_by_type(df2, "individual_LDstats_mean_nonzero_by_type.png")

ttest_IndvLD_by_type(df2)

plot_IndvLD_by_comparison(df2, "individual_LDstats_mean_nonzero_by_type.png")

plot_IndvLD_by_comparison_experimental(df2, "individual_LDstats_mean_nonzero_sorted.png")

plot_IndvLD_ecdf_by_comparison(df2, "individual_LDstats_mean_nonzero_ecdf.png")

