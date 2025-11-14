import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_branch_stats(csvfile, title, stat1, stat2, colors, inversion):
    """
    Plots the longest branch length and normalized longest branch length from the csv output of windowed_hclustering.py
    If there is a known inversion in this region, it can be highlighted by setting the inversion argument to a tuple of 
    2 breakpoints. If there is no inversion set to None.

    Arguments:
    csvfile (str): path to datafile to plot from
    title (str): title of the plot and name of the png file saved by the function
    stat1 (str): name of column with data for the first subplot and will also be used for the y axis label
    stat2 (str): name of column with data for the second subplot and will also be used for the y axis label
    colors (tuple of 2 str): tuple of color names for plotting
    inversion (tuple of 2 int or Nonetype): tuple of inversion breakpoints. If None, no inversion breakpoints will be plotted
    """

    #opening file
    df = pd.read_csv(csvfile)

    #creating plot
    fig, axs = plt.subplots(2, 1, figsize=(15, 5))

    #plot dip stat
    axs[0].plot(df["POS"], df[stat1], color=colors[0])
    if inversion:
        axs[0].axvspan(inversion[0], inversion[1], color='yellow', alpha=0.4)
    #axs[0].set_title("Normalized Longest Branch Length")
    axs[0].set_ylabel(stat1.replace("_", " "))
    axs[0].ticklabel_format(style='plain')
    axs[0].set_ylim(0.25,1)

    #plot pve
    axs[1].plot(df["POS"], df[stat2], color=colors[1])
    if inversion:
        axs[1].axvspan(inversion[0], inversion[1], color='yellow', alpha=0.4)
    #axs[1].set_title("PC1 Percent Variance Explained")
    axs[1].set_ylabel(stat2.replace("_", " "))
    axs[1].set_xlabel("Position")
    axs[1].ticklabel_format(style='plain')
    axs[1].set_ylim(0,0.5)

    plt.suptitle(title)

    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{title}.png")
    plt.clf()







inv_map = {
    "chr9": (113105263, 113119310),
    "chr19": (21647331, 22062459),
    "chr17": (36357258, 37957978),
    "chr15.1": (30077909, 32607507),
    "chr15.2": (22770522, 28852548)
}


# file_list = [
#     "top_branch_lengths200_chr15.1_1MB_inv_1MB_4529599bp.csv", "top_branch_lengths200_chr15.1_control1_4529599bp.csv", "top_branch_lengths200_chr15.1_control2_4529599bp.csv",
#     "top_branch_lengths200_chr15.2_1MB_inv_1MB_8082027bp.csv", "top_branch_lengths200_chr15.2_control1_8082027bp.csv", "top_branch_lengths200_chr15.2_control2_8082027bp.csv",
#     "top_branch_lengths200_chr17_1MB_inv_1MB_3600721bp.csv", "top_branch_lengths200_chr17_control1_3600721bp.csv", "top_branch_lengths200_chr17_control2_3600721bp.csv",
#     "top_branch_lengths200_chr19_1MB_inv_1MB_2415129bp.csv", "top_branch_lengths200_chr19_control1_2415129bp.csv", "top_branch_lengths200_chr19_control2_2415129bp.csv",
#     "top_branch_lengths200_chr9_1MB_inv_1MB_2014048bp.csv", "top_branch_lengths200_chr9_control1_2014048bp.csv", "top_branch_lengths200_chr9_control2_2014048bp.csv",
#     "top_branch_lengths400_chr15.1_1MB_inv_1MB_4529599bp.csv", "top_branch_lengths400_chr15.1_control1_4529599bp.csv", "top_branch_lengths400_chr15.1_control2_4529599bp.csv",
#     "top_branch_lengths400_chr15.2_1MB_inv_1MB_8082027bp.csv", "top_branch_lengths400_chr15.2_control1_8082027bp.csv", "top_branch_lengths400_chr15.2_control2_8082027bp.csv",
#     "top_branch_lengths400_chr17_1MB_inv_1MB_3600721bp.csv", "top_branch_lengths400_chr17_control1_3600721bp.csv", "top_branch_lengths400_chr17_control2_3600721bp.csv",
#     "top_branch_lengths400_chr19_1MB_inv_1MB_2415129bp.csv", "top_branch_lengths400_chr19_control1_2415129bp.csv", "top_branch_lengths400_chr19_control2_2415129bp.csv",
#     "top_branch_lengths400_chr9_1MB_inv_1MB_2014048bp.csv", "top_branch_lengths400_chr9_control1_2014048bp.csv", "top_branch_lengths400_chr9_control2_2014048bp.csv",
# ]


file_list = [
    "LongestAvg_branch_lengths_200SNP_windows_chr15.1_inv.csv", "LongestAvg_branch_lengths_200SNP_windows_chr15.1_control1.csv", "LongestAvg_branch_lengths_200SNP_windows_chr15.1_control2.csv",
    "LongestAvg_branch_lengths_200SNP_windows_chr15.2_inv.csv", "LongestAvg_branch_lengths_200SNP_windows_chr15.2_control1.csv", "LongestAvg_branch_lengths_200SNP_windows_chr15.2_control2.csv",
    "LongestAvg_branch_lengths_200SNP_windows_chr17_inv.csv", "LongestAvg_branch_lengths_200SNP_windows_chr17_control1.csv", "LongestAvg_branch_lengths_200SNP_windows_chr17_control2.csv",
    "LongestAvg_branch_lengths_200SNP_windows_chr19_inv.csv", "LongestAvg_branch_lengths_200SNP_windows_chr19_control1.csv", "LongestAvg_branch_lengths_200SNP_windows_chr19_control2.csv",
    "LongestAvg_branch_lengths_200SNP_windows_chr9_inv.csv", "LongestAvg_branch_lengths_200SNP_windows_chr9_control1.csv", "LongestAvg_branch_lengths_200SNP_windows_chr9_control2.csv",
    "LongestAvg_branch_lengths_400SNP_windows_chr15.1_inv.csv", "LongestAvg_branch_lengths_400SNP_windows_chr15.1_control1.csv", "LongestAvg_branch_lengths_400SNP_windows_chr15.1_control2.csv",
    "LongestAvg_branch_lengths_400SNP_windows_chr15.2_inv.csv", "LongestAvg_branch_lengths_400SNP_windows_chr15.2_control1.csv", "LongestAvg_branch_lengths_400SNP_windows_chr15.2_control2.csv",
    "LongestAvg_branch_lengths_400SNP_windows_chr17_inv.csv", "LongestAvg_branch_lengths_400SNP_windows_chr17_control1.csv", "LongestAvg_branch_lengths_400SNP_windows_chr17_control2.csv",
    "LongestAvg_branch_lengths_400SNP_windows_chr19_inv.csv", "LongestAvg_branch_lengths_400SNP_windows_chr19_control1.csv", "LongestAvg_branch_lengths_400SNP_windows_chr19_control2.csv",
    "LongestAvg_branch_lengths_400SNP_windows_chr9_inv.csv", "LongestAvg_branch_lengths_400SNP_windows_chr9_control1.csv", "LongestAvg_branch_lengths_400SNP_windows_chr9_control2.csv",
]


for file in file_list:
    chrom = file.split("_")[-2]
    window_size = file.split("_")[3].replace("SNP", "")
    if "inv" in file:
        inverted = "inv"
        breakpoints = inv_map[chrom]
    else:
        inverted = file.split("_")[-1].replace(".csv", "")
        breakpoints = None
    plottitle = f"{chrom} {inverted} in {window_size} SNP windows"


    plot_branch_stats(file, plottitle, "Longest_branch_length", "Average_branch_length", ("green", "blue"), breakpoints)