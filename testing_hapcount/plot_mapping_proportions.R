library(ggplot2)
library(forcats)
library(patchwork)


setwd("~/Desktop/Phadnis_Lab/human_inversions/testing_hapcount_classification_power/counting_intersects")


calculate_mapping_proportion <- function(bedfile){
    ###Function takes in a bedfile if intersect count and outputs the proportion of spans that mapped to the intersecting file
    
    #opening bed
    df = read.table(bedfile, sep="\t", col.names=c("CHROM", "START", "END", "COUNT"))
    
    #did the query span map or not
    df$mapped <- as.integer(df$COUNT > 0)
    
    #total
    tot <- nrow(df)
    
    #number mapped
    num_mapped <- sum(df$mapped)
    
    #proportion
    prop <- num_mapped / tot
    
    return(prop)
    
}



###Plotting proportions of SD and Invs Mapped to each other
prop_SDs_mapped_to_Invs <- calculate_mapping_proportion("SDs_count_mapping_Invs.bed")
prop_Invs_mapped_to_SDs <- calculate_mapping_proportion("Invs_count_mapping_SDs.bed")

res = 300
scale = 2.5
png("SD_Inv_proprtions.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
par(mfrow = c(1, 2))
pie(c(prop_SDs_mapped_to_Invs, 1-prop_SDs_mapped_to_Invs), labels=c("Inversions", "No Inversions"),
    col=c("plum", "lightskyblue"), main=paste0(round(prop_SDs_mapped_to_Invs*100, 1), "% of SDs\nmapped to Inversions"))
pie(c(prop_Invs_mapped_to_SDs, 1-prop_Invs_mapped_to_SDs), labels=c("SDs", "No SDs"),
    col=c("plum", "lightskyblue"), main=paste0(round(prop_Invs_mapped_to_SDs*100, 1), "% of Inversions\nmapped to SDs"))
dev.off()
#resetting plot
par(mfrow = c(1, 1))



###Plotting proportions of SD and InvBrks Mapped to each other
prop_SDs_mapped_to_InvBrks <- calculate_mapping_proportion("SDs_count_mapping_InvBrks.bed")
prop_InvBrks_mapped_to_SDs <- calculate_mapping_proportion("InvBrks_count_mapping_SDs.bed")

res = 300
scale = 2.5
png("SD_InvBrks_proprtions.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
par(mfrow = c(1, 2))
pie(c(prop_SDs_mapped_to_InvBrks, 1-prop_SDs_mapped_to_InvBrks), labels=c("Inv\nBrkpnts", "No Inv\nBrkpnts"),
    col=c("plum", "lightskyblue"), main=paste0(round(prop_SDs_mapped_to_InvBrks*100, 1), "% of SDs\nmapped to Inv Breakpoints"))
pie(c(prop_InvBrks_mapped_to_SDs, 1-prop_InvBrks_mapped_to_SDs), labels=c("SDs", "No SDs"),
    col=c("plum", "lightskyblue"), main=paste0(round(prop_InvBrks_mapped_to_SDs*100, 1), "% of Inv Breakpoints\nmapped to SDs"))
dev.off()
#resetting plot
par(mfrow = c(1, 1))



###Plotting proportion of Peaks mapped to SDs and Inversions

#making dfs

f1 <- "Peaks_lowess90_100Kbwin_count_mapping_SDs.bed"
f2 <- "Peaks_lowess95_100Kbwin_count_mapping_SDs.bed"
f3 <- "Peaks_findpeaks90_100Kbwin_count_mapping_SDs.bed"
f4 <- "Peaks_findpeaks95_100Kbwin_count_mapping_SDs.bed"
f5 <- "Peaks_lowessfindpeaks90_100Kbwin_count_mapping_SDs.bed"
f6 <- "Peaks_lowess90_10Kbwin_count_mapping_SDs.bed"
f7 <- "Peaks_lowess95_10Kbwin_count_mapping_SDs.bed"
f8 <- "Peaks_findpeaks90_10Kbwin_count_mapping_SDs.bed"
f9 <- "Peaks_findpeaks95_10Kbwin_count_mapping_SDs.bed"
f10 <- "Peaks_lowessfindpeaks90_10Kbwin_count_mapping_SDs.bed"

peak_SD <- data.frame(
    Method = c(
        "P90 Lowess\n100Kb windows",
        "P95 Lowess\n100Kb windows",
        "P90 find_peaks\n100Kb windows",
        "P95 find_peaks\n100Kb windows",
        "P90 Lowess & find_peaks\n100Kb windows",
        "P90 Lowess\n10Kb windows",
        "P95 Lowess\n10Kb windows",
        "P90 find_peaks\n10Kb windows",
        "P95 find_peaks\n10Kb windows",
        "P90 Lowess & find_peaks\n10Kb windows"
    ),
    Proportion = c(
        calculate_mapping_proportion(f1),
        calculate_mapping_proportion(f2),
        calculate_mapping_proportion(f3),
        calculate_mapping_proportion(f4),
        calculate_mapping_proportion(f5),
        calculate_mapping_proportion(f6),
        calculate_mapping_proportion(f7),
        calculate_mapping_proportion(f8),
        calculate_mapping_proportion(f9),
        calculate_mapping_proportion(f10)
    )
)



f11 <- "Peaks_lowess90_100Kbwin_count_mapping_Invs.bed"
f12 <- "Peaks_lowess95_100Kbwin_count_mapping_Invs.bed"
f13 <- "Peaks_findpeaks90_100Kbwin_count_mapping_Invs.bed"
f14 <- "Peaks_findpeaks95_100Kbwin_count_mapping_Invs.bed"
f15 <- "Peaks_lowessfindpeaks90_100Kbwin_count_mapping_Invs.bed"
f16 <- "Peaks_lowess90_10Kbwin_count_mapping_Invs.bed"
f17 <- "Peaks_lowess95_10Kbwin_count_mapping_Invs.bed"
f18 <- "Peaks_findpeaks90_10Kbwin_count_mapping_Invs.bed"
f19 <- "Peaks_findpeaks95_10Kbwin_count_mapping_Invs.bed"
f20 <- "Peaks_lowessfindpeaks90_10Kbwin_count_mapping_Invs.bed"


peak_inv <- data.frame(
    Method = c(
        "P90 Lowess\n100Kb windows",
        "P95 Lowess\n100Kb windows",
        "P90 find_peaks\n100Kb windows",
        "P95 find_peaks\n100Kb windows",
        "P90 Lowess & find_peaks\n100Kb windows",
        "P90 Lowess\n10Kb windows",
        "P95 Lowess\n10Kb windows",
        "P90 find_peaks\n10Kb windows",
        "P95 find_peaks\n10Kb windows",
        "P90 Lowess & find_peaks\n10Kb windows"
    ),
    Proportion = c(
        calculate_mapping_proportion(f11),
        calculate_mapping_proportion(f12),
        calculate_mapping_proportion(f13),
        calculate_mapping_proportion(f14),
        calculate_mapping_proportion(f15),
        calculate_mapping_proportion(f16),
        calculate_mapping_proportion(f17),
        calculate_mapping_proportion(f18),
        calculate_mapping_proportion(f19),
        calculate_mapping_proportion(f20)
    )
)




p1 <- ggplot(data=peak_SD, aes(x=fct_reorder(Method, -Proportion), y=Proportion, fill=Method)) +
    geom_col() +
    ggtitle("Proportion of Peaks Mapped to SDs") + xlab("Peak Finding Method") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())



p2 <- ggplot(data=peak_inv, aes(x=fct_reorder(Method, -Proportion), y=Proportion, fill=Method)) +
    geom_col() +
    ggtitle("Proportion of Peaks Mapped to Inversions") + xlab("Peak Finding Method") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
    

res = 300
scale = 3.5
png("prop_peaks_mapped_to_SDs_Invs.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
p1 + p2
dev.off()












###Plotting proportion of Peaks mapped to SDs and Inversion BREAKPOINTS

#making dfs

f1 <- "Peaks_lowess90_100Kbwin_count_mapping_SDs.bed"
f2 <- "Peaks_lowess95_100Kbwin_count_mapping_SDs.bed"
f3 <- "Peaks_findpeaks90_100Kbwin_count_mapping_SDs.bed"
f4 <- "Peaks_findpeaks95_100Kbwin_count_mapping_SDs.bed"
f5 <- "Peaks_lowessfindpeaks90_100Kbwin_count_mapping_SDs.bed"
f6 <- "Peaks_lowess90_10Kbwin_count_mapping_SDs.bed"
f7 <- "Peaks_lowess95_10Kbwin_count_mapping_SDs.bed"
f8 <- "Peaks_findpeaks90_10Kbwin_count_mapping_SDs.bed"
f9 <- "Peaks_findpeaks95_10Kbwin_count_mapping_SDs.bed"
f10 <- "Peaks_lowessfindpeaks90_10Kbwin_count_mapping_SDs.bed"

peak_SD <- data.frame(
    Method = c(
        "P90 Lowess\n100Kb windows",
        "P95 Lowess\n100Kb windows",
        "P90 find_peaks\n100Kb windows",
        "P95 find_peaks\n100Kb windows",
        "P90 Lowess & find_peaks\n100Kb windows",
        "P90 Lowess\n10Kb windows",
        "P95 Lowess\n10Kb windows",
        "P90 find_peaks\n10Kb windows",
        "P95 find_peaks\n10Kb windows",
        "P90 Lowess & find_peaks\n10Kb windows"
    ),
    Proportion = c(
        calculate_mapping_proportion(f1),
        calculate_mapping_proportion(f2),
        calculate_mapping_proportion(f3),
        calculate_mapping_proportion(f4),
        calculate_mapping_proportion(f5),
        calculate_mapping_proportion(f6),
        calculate_mapping_proportion(f7),
        calculate_mapping_proportion(f8),
        calculate_mapping_proportion(f9),
        calculate_mapping_proportion(f10)
    )
)



f11 <- "Peaks_lowess90_100Kbwin_count_mapping_InvBrks.bed"
f12 <- "Peaks_lowess95_100Kbwin_count_mapping_InvBrks.bed"
f13 <- "Peaks_findpeaks90_100Kbwin_count_mapping_InvBrks.bed"
f14 <- "Peaks_findpeaks95_100Kbwin_count_mapping_InvBrks.bed"
f15 <- "Peaks_lowessfindpeaks90_100Kbwin_count_mapping_InvBrks.bed"
f16 <- "Peaks_lowess90_10Kbwin_count_mapping_InvBrks.bed"
f17 <- "Peaks_lowess95_10Kbwin_count_mapping_InvBrks.bed"
f18 <- "Peaks_findpeaks90_10Kbwin_count_mapping_InvBrks.bed"
f19 <- "Peaks_findpeaks95_10Kbwin_count_mapping_InvBrks.bed"
f20 <- "Peaks_lowessfindpeaks90_10Kbwin_count_mapping_InvBrks.bed"


peak_inv <- data.frame(
    Method = c(
        "P90 Lowess\n100Kb windows",
        "P95 Lowess\n100Kb windows",
        "P90 find_peaks\n100Kb windows",
        "P95 find_peaks\n100Kb windows",
        "P90 Lowess & find_peaks\n100Kb windows",
        "P90 Lowess\n10Kb windows",
        "P95 Lowess\n10Kb windows",
        "P90 find_peaks\n10Kb windows",
        "P95 find_peaks\n10Kb windows",
        "P90 Lowess & find_peaks\n10Kb windows"
    ),
    Proportion = c(
        calculate_mapping_proportion(f11),
        calculate_mapping_proportion(f12),
        calculate_mapping_proportion(f13),
        calculate_mapping_proportion(f14),
        calculate_mapping_proportion(f15),
        calculate_mapping_proportion(f16),
        calculate_mapping_proportion(f17),
        calculate_mapping_proportion(f18),
        calculate_mapping_proportion(f19),
        calculate_mapping_proportion(f20)
    )
)




p1 <- ggplot(data=peak_SD, aes(x=fct_reorder(Method, -Proportion), y=Proportion, fill=Method)) +
    geom_col() +
    ggtitle("Proportion of Peaks Mapped to SDs") + xlab("Peak Finding Method") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())



p2 <- ggplot(data=peak_inv, aes(x=fct_reorder(Method, -Proportion), y=Proportion, fill=Method)) +
    geom_col() +
    ggtitle("Proportion of Peaks Mapped to Inv Breakpoints") + xlab("Peak Finding Method") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())


res = 300
scale = 3.5
png("prop_peaks_mapped_to_SDs_InvBrks.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
p1 + p2
dev.off()





















#better plot for viewing



###Plotting proportion of Peaks mapped to SDs and Inversion BREAKPOINTS

#making dfs

f1 <- "Peaks_lowess90_100Kbwin_count_mapping_SDs.bed"
f2 <- "Peaks_lowess95_100Kbwin_count_mapping_SDs.bed"
f3 <- "Peaks_findpeaks90_100Kbwin_count_mapping_SDs.bed"
f4 <- "Peaks_findpeaks95_100Kbwin_count_mapping_SDs.bed"
f5 <- "Peaks_lowessfindpeaks90_100Kbwin_count_mapping_SDs.bed"


peak_SD <- data.frame(
    Method = c(
        "P90 Lowess",
        "P95 Lowess",
        "P90 find_peaks",
        "P95 find_peaks",
        "P90 Lowess & find_peaks"
    ),
    Proportion = c(
        calculate_mapping_proportion(f1),
        calculate_mapping_proportion(f2),
        calculate_mapping_proportion(f3),
        calculate_mapping_proportion(f4),
        calculate_mapping_proportion(f5)
    )
)



f11 <- "Peaks_lowess90_100Kbwin_count_mapping_InvBrks.bed"
f12 <- "Peaks_lowess95_100Kbwin_count_mapping_InvBrks.bed"
f13 <- "Peaks_findpeaks90_100Kbwin_count_mapping_InvBrks.bed"
f14 <- "Peaks_findpeaks95_100Kbwin_count_mapping_InvBrks.bed"
f15 <- "Peaks_lowessfindpeaks90_100Kbwin_count_mapping_InvBrks.bed"


peak_inv <- data.frame(
    Method = c(
        "P90 Lowess",
        "P95 Lowess",
        "P90 find_peaks",
        "P95 find_peaks",
        "P90 Lowess & find_peaks"
    ),
    Proportion = c(
        calculate_mapping_proportion(f11),
        calculate_mapping_proportion(f12),
        calculate_mapping_proportion(f13),
        calculate_mapping_proportion(f14),
        calculate_mapping_proportion(f15)
    )
)




p1 <- ggplot(data=peak_SD, aes(x=fct_reorder(Method, -Proportion), y=Proportion, fill=Method)) +
    geom_col() +
    ggtitle("Proportion of Peaks Mapped to SDs") + xlab("Peak Finding Method") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())



p2 <- ggplot(data=peak_inv, aes(x=fct_reorder(Method, -Proportion), y=Proportion, fill=Method)) +
    geom_col() +
    ggtitle("Proportion of Peaks Mapped to Inv Breakpoints") + xlab("Peak Finding Method") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())


res = 300
scale = 3
png("BETTER_PLOT_prop_peaks_mapped_to_SDs_InvBrks.png", height = 3 * scale * res, width = 4 * scale * res, res = res)
p1 + p2
dev.off()
