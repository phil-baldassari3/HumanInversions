library(ggplot2)
library(RColorBrewer)



setwd("~/Desktop/Phadnis_Lab/human_inversions/testing_hapcount_classification_power/bootstrap_LD")



plot_LD_between_2regions <- function(plink_output, region1_label, region2_label){
    ### Function takes a plink output from my python program `bootstrap_LD.py`
    ### and plots |D'| between region 1 and region 2
    
    #opening data
    df <- read.csv(plink_output)
    
    #plotting heatmap
    #plotting
    outpng <- paste0(region1_label, "__TO__", region2_label, ".png")
    res = 300
    scale = 1
    png(outpng, height = 3 * scale * res, width = 4 * scale * res, res = res)
    p <- ggplot(data=df, aes(x=idx_A, y=idx_B, fill=Dprime)) +
        geom_tile() + scale_fill_distiller(palette = "YlOrRd", direction=1, name="|D'|") +
        xlab(paste(region1_label, "SNP Index")) + ylab(paste(region2_label, "SNP Index")) +
        theme_classic()
    print(p)
    dev.off()

    
}





plot_LD_between_2regions(
    "peak49.chr21_26310002_26620000.invbrk0.NA_NA_NA__TO__peak50.chr21_27760002_28180000.invbrk0.NA_NA_NA.ld", 
    "Chr21_peak49_None", "Chr21_peak50_None"
)

plot_LD_between_2regions(
    "peak40.chr17_16380002_17090000.invbrk9.chr17_16823491_NA__TO__peak41.chr17_18270002_19250000.invbrk9.chr17_NA_18384190.ld", 
    "Chr17_peak40_INV", "Chr17_peak41_INV"
)

plot_LD_between_2regions(
    "peak39.chr17_15310002_15950000.invbrk0.NA_NA_NA__TO__peak40.chr17_16380002_17090000.invbrk9.chr17_16823491_NA.ld", 
    "Chr17_peak39_None", "Chr17_peak40_INV"
)

plot_LD_between_2regions(
    "peak38.chr17_13680002_14470000.invbrk0.NA_NA_NA__TO__peak41.chr17_18270002_19250000.invbrk9.chr17_NA_18384190.ld", 
    "Chr17_peak38_None", "Chr17_peak41_INV"
)

plot_LD_between_2regions(
    "peak38.chr17_13680002_14470000.invbrk0.NA_NA_NA__TO__peak40.chr17_16380002_17090000.invbrk9.chr17_16823491_NA.ld", 
    "Chr17_peak38_None", "Chr17_peak40_INV"
)

plot_LD_between_2regions(
    "peak38.chr17_13680002_14470000.invbrk0.NA_NA_NA__TO__peak39.chr17_15310002_15950000.invbrk0.NA_NA_NA.ld", 
    "Chr17_peak38_None", "Chr17_peak39_None"
)

plot_LD_between_2regions(
    "peak31.chr15_31950002_32820000.invbrk6.chr15_NA_32153204__TO__peak32.chr15_34250002_34800000.invbrk0.NA_NA_NA.ld", 
    "Chr15_peak31_INV", "Chr15_peak32_None"
)

plot_LD_between_2regions(
    "peak30.chr15_29990002_30910000.invbrk6.chr15_30618104_NA__TO__peak31.chr15_31950002_32820000.invbrk6.chr15_NA_32153204.ld", 
    "Chr15_peak30_INV", "Chr15_peak31_INV"
)

plot_LD_between_2regions(
    "peak29.chr15_28226002_28784000.invbrk5.chr15_NA_28389868__TO__peak30.chr15_29990002_30910000.invbrk6.chr15_30618104_NA.ld", 
    "Chr15_peak29_INV5", "Chr15_peak30_INV6"
)

plot_LD_between_2regions(
    "peak28.chr15_23010002_23361000.invbrk5.chr15_23345460_NA__TO__peak29.chr15_28226002_28784000.invbrk5.chr15_NA_28389868.ld", 
    "Chr15_peak28_INV", "Chr15_peak29_INV"
)

plot_LD_between_2regions(
    "peak26.chr13_36490002_36980000.invbrk0.NA_NA_NA__TO__peak27.chr13_45130002_47250000.invbrk0.NA_NA_NA.ld", 
    "Chr13_peak26_None", "Chr13_peak27_None"
)

plot_LD_between_2regions(
    "peak25.chr13_35560002_36150000.invbrk0.NA_NA_NA__TO__peak26.chr13_36490002_36980000.invbrk0.NA_NA_NA.ld", 
    "Chr13_peak25_None", "Chr13_peak26_None"
)

plot_LD_between_2regions(
    "peak19.chr8_11760002_12860000.invbrk3.chr8_NA_12598379__TO__peak20.chr8_16610002_17900000.invbrk0.NA_NA_NA.ld", 
    "Chr_peak19_INV", "Chr_peak20_None"
)

plot_LD_between_2regions(
    "peak18.chr8_6790002_8380000.invbrk3.chr8_7301025_NA__TO__peak20.chr8_16610002_17900000.invbrk0.NA_NA_NA.ld",
    "Chr8_peak18_INV", "Chr8_peak20_None"
)

plot_LD_between_2regions(
    "peak18.chr8_6790002_8380000.invbrk3.chr8_7301025_NA__TO__peak19.chr8_11760002_12860000.invbrk3.chr8_NA_12598379.ld",
    "Chr8_peak18_INV", "Chr8_peak19_INV"
)

plot_LD_between_2regions(
    "peak17.chr8_3560002_4560000.invbrk0.NA_NA_NA__TO__peak19.chr8_11760002_12860000.invbrk3.chr8_NA_12598379.ld",
    "Chr8_peak17_None", "Chr8_peak19_INV"
)

plot_LD_between_2regions(
    "peak17.chr8_3560002_4560000.invbrk0.NA_NA_NA__TO__peak18.chr8_6790002_8380000.invbrk3.chr8_7301025_NA.ld",
    "Chr8_peak17_None", "Chr8_peak18_INV"
)


plot_LD_between_2regions(
    "peak7.chr2_111080002_112210000.invbrk2.chr2_NA_111255403__TO__peak8.chr2_120190002_121520000.invbrk0.NA_NA_NA.ld", 
    "Chr2_peak7_INV", "Chr2_peak8_None"
)

plot_LD_between_2regions(
    "peak6.chr2_86590002_88340000.invbrk2.chr2_87987172_NA__TO__peak7.chr2_111080002_112210000.invbrk2.chr2_NA_111255403.ld", 
    "Chr2_peak6_INV", "Chr2_peak7_INV"
)

plot_LD_between_2regions(
    "peak5.chr2_76140002_78120000.invbrk0.NA_NA_NA__TO__peak6.chr2_86590002_88340000.invbrk2.chr2_87987172_NA.ld", 
    "Chr2_peak5_None", "Chr2_peak6_INV"
)
