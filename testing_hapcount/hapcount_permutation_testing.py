"""
This is version 3 of this script. It only tests for a significant difference in 
proportion of unique haplotypes between coordinate windows and whole genome.
There is also an additional function that allows you to convert a bedfile of 
coordinates to a list of tuples to be input into the main function.

Note that averaging was done as a whole instead of computing an average of averages! You should probably fix this later!
#I implemented a (temporary) fix for this noted with "###"
"""


import os,sys
import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev
import scipy.stats
import random
import matplotlib.pyplot as plt


random.seed(349)



def _ztest(top_value, null_set):
    """
    Performs a two-tailed Z test between the set of interest and the null set

    :param top_value: integer or float for the average value from the set of interest
    :type top_value: float
    :param null_set: list of integers or floats that of average values in the null set from the permutation
    :type null_set: list of floats
    
    :returns: z-score and p-value
    :rtype: float, float
    """

    #to avoid division by zero error but this means that the z test is meaningless
    if stdev(null_set) == 0:
        zscore = 0.00
        pvalue = 1.00

    else:   
        zscore = round((top_value - mean(null_set)) / stdev(null_set), 2)
        pvalue = round(2 * scipy.stats.norm.sf(zscore), 2)

    return zscore, pvalue



def _avg_from_windows_within_coordinates(df, coordinates, val_param):
    """
    Filters a given dataframe by the specified genomic coordinates and outputs the average of specified values from the dataframe

    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param coordinates: list of tuples with regions to include [("chrom", start, end)...]
    :type coordinates: list of tuples of str, int, and int
    :param val_param: column name for which to grab values after filtering
    :type val_param: str

    :returns: average of values after filtering and a list the number of windows included in each breakpoint span
    :rtype: float, list
    """

    val_list = []
    num_win_per_brk = []

    #bool masking rows for windows entirely or partially overlapping with regions to keep
    for chrom, region_start, region_end in coordinates:
        keep_rows = ((df["CHROM"] == chrom) & (df["START"] < region_end) & (df["END"] > region_start))
        
        #filter data
        filtered_df = df.loc[keep_rows]

        #grabbing values
        brk_val_list = filtered_df[val_param].to_list()

        ###
        #check if there is any data at this coordinate, if not just skip it
        if len(brk_val_list) == 0:
            continue

        #grabbing number of windows in breakpoint
        num_win_per_brk.append(len(brk_val_list))

        #add brk vals to big vals list
        #val_list += brk_val_list
        ####
        val_list.append(mean(brk_val_list))

    #averaging values
    avg_val = mean(val_list)

    return avg_val, num_win_per_brk



def _sample_windows_for_1_brkpnt(df, num_windows, data_param):
    """
    Helper function that randomly samples a genomic span of a target length. The function makes sure that the span is within the length of the chromosome.
    
    :param df: windowed genomic data
    :param num_windows: target number of windows in a row to sample from the windowed data
    :param data_param: column name for values to sample and average

    :type data_corrdinates: [(str, int, int)]
    :type length_of_span: int
    :type data_param: str

    :returns: list of data from sampled windows
    :rtype: list
    """

    sampled_data = []

    while len(sampled_data) == 0:
        #sample window start coordinate and subset dataframe for range of windows
        sampled_idx = random.randrange(len(df))
        sampled_df = df.iloc[sampled_idx:sampled_idx+num_windows]
        #check if target end point is within the chromosome, if not, sample again
        if sampled_df["CHROM"].iloc[0] == sampled_df["CHROM"].iloc[-1]:
            sampled_data = sampled_df[data_param].to_list()
            break

    return sampled_data


def sample_windows_for_brkpnt_set(data, window_counts, colname):
    """
    Function samples random genomic spans corresponding to each given breakpoint region for a given number of permutations using `_randomly_sample_span`.
    
    :param data: windowed genomic data
    :param window_counts: list of number of windows in each breakpoint found in `_avg_from_windows_within_coordinates` and used for sampling
    :param colname: column name for values to sample and average

    :type data: DataFrame
    :type breakpoints: [(str, int, int)]
    :type colname: str

    :returns: randomly sampled data corresponding to the breakpoint set
    :rtype: list
    """

    #empty list that will hold sampled data corresponding to the whole breakpoint set
    sampled_data_set = []

    #for each breakpoint sample corresponding random span
    for count in window_counts:

        if count == 0:
            continue

        #grab sampled data for corresponding breakpoint
        sampled_vals = _sample_windows_for_1_brkpnt(data, count, colname)

        #sampled_data_set += sampled_vals
        ###
        sampled_data_set.append(mean(sampled_vals))


    return sampled_data_set



def _permute_data(df, param_to_permute, list_of_num_windows, permutations):
    """
    Permutes data a given number of times, each time grabbing a random sampled windows corresponding to each  
    breakpoint length, finding the average of all these windows, and appending to a list to generate a null distribution.
    
    :param df: dataframe of windowed proportion of unique haplotypes
    :param param_to_permute: column name for values to sample and average
    :param list_of_num_windows: list of number of windows in each breakpoint found in `_avg_from_windows_within_coordinates` and used for sampling
    :param permutations: number of random permutations

    :returns: lists of averaged values, one for each permutation thus generating a null distribution
    :rtype: list of floats
    """

    #empty list for null distribution
    null_dist = []

    #permutating
    for i in range(permutations):
        if i % 1000 == 0:
            print(f"Permutation #{i}")

        #sampling random coordinates cooresponding to each breakpoint
        sampled_data_vals = sample_windows_for_brkpnt_set(df, list_of_num_windows, param_to_permute)
        avg_of_sample = mean(sampled_data_vals)
        null_dist.append(avg_of_sample)


    return null_dist



def _plot_z_test(null_list, val_of_interest, param_name, plot_title, zscore, pvalue):
    """
    Plots the null distribution in a histogram and the value of interest as a vertical line. Saves a png of plot
    
    :param null_list: null distribution from _permute_data
    :type null_list: list
    :param val_of_interest: value to be compared to the null distribution
    :type val_of_interest: float
    :param param_name: name of parameter being plotted, used for x-axis label
    :type param_name: str
    :param plot_title: title of plot, also used as the filename for output png
    :type plot_title: str
    :param zscore: z-score from _ztest, used to display as text on plot
    :type zscore: float
    :param pvalue: p-value from _ztest, used to display as text on plot
    :type pvalue: float
    """

    #number of bins using Freedman–Diaconis rule
    q25, q75 = np.percentile(null_list, [25, 75])
    bin_width = 2 * (q75 - q25) * len(null_list) ** (-1/3)
    
    if bin_width == 0:
        bins = 10
    else:
        try:
            bins = round((max(null_list) - min(null_list)) / bin_width)
        except OverflowError:
            bins = 10

    #plotting
    plt.figure(figsize=(8,6))
    plt.hist(null_list, density=True, bins=bins)
    plt.ylabel('Density')
    plt.xlabel(param_name.replace("_", " "))
    plt.title(plot_title)
  
    plt.axvline(x=val_of_interest, color='r', lw=3)

    plt.text(0.5, 0.8, f"z-score = {round(zscore, 3)}\np-value = {pvalue}", 
             horizontalalignment='center', verticalalignment='bottom', bbox=dict(facecolor='yellow', alpha=0.25),
             transform=plt.gca().transAxes)

    plt.savefig(f"{plot_title}.png", bbox_inches='tight')
    #plt.show()
    plt.clf()



def breakpoint_coordinates_from_bed(bedfile, bp):
    """
    Function takes in a bedfile of coordinates and outputs a list of breakpoint coordinates.
    
    :param bedfile: .bed file of inversion spans
    :param bp: number of basepairs to include on either side of the breakpoint in the output coordinates
    :param makeuniq: optional argument to remove duplicate spans in output

    :type bedfile: str
    :type bp: int
    :type makeuniq: bool=False

    :returns: list of tuples each with CHROM, START, and END defining breakpoint coordinates
    :rtype: [(str, int, int)]
    """

    print("Reading in inversion coordinates...")

    #output
    coords = []

    #open bedfile and iterate through lines
    with open(bedfile, "r") as bed:
        for line in bed:
            linels = line.strip().split("\t")

            #grab coords
            chrom = linels[0]
            start1 = int(linels[1]) - bp
            end1 = int(linels[1]) + bp
            start2 = int(linels[2]) - bp
            end2 = int(linels[2]) + bp

            #package into tuples
            tup1 = (chrom, start1, end1)
            tup2 = (chrom, start2, end2)

            #append to coords
            coords.append(tup1)
            coords.append(tup2)


    return coords



def SD_coordinates_from_bed(bedfile, bp):
    """
    Function takes in a bedfile of SDs and outputs a list of coordinates.
    
    :param bedfile: .bed file of spans
    :param bp: number of basepairs to include on either side of the breakpoint in the output coordinates

    :type bedfile: str
    :type bp: int

    :returns: list of tuples each with CHROM, START, and END defining breakpoint coordinates
    :rtype: [(str, int, int)]
    """

    print("Reading in SD coordinates...")

    #output
    coords = []

    #open bedfile and iterate through lines
    with open(bedfile, "r") as bed:
        for line in bed:
            linels = line.strip().split("\t")

            #grab coords
            chrom = linels[0]
            start = int(linels[1]) - bp
            end = int(linels[2]) + bp

            #package into tuple
            tup = (chrom, start, end)

            #append to coords
            coords.append(tup)

    #make unique
    coords = list(set(coords))

    return coords




def run_ztest_on_coordinates(bedgraph, breakpoint_coordinates, num_permutations, plottitle):
    """
    Main function to run permutation Z-test on last bedgraph column (Proportion of unique haplotypes) to test for a significant difference
    between its values in breakpoints versus genome wide
    
    :param bedgraph: bedgraph file of windowed data. There shouldnt be a header line but the columns should be CHROM, START, END, Proportion of Unique Haplotypes
    :param coordinates_to_extract: list of tuples with regions to include [("chrom", start, end)...]
    :param num_permutations: number of permutations
    :param plottitle: Title of plot and name of filename

    :type bedgraph: str
    :type coordinates_to_extract: list of tuples of str, int, and int
    :type num_permutations: int
    :type plottitle: str
    """

    print("Running Permutation...")

    data = pd.read_csv(bedgraph, sep="\t", header=None, names=["CHROM", "START", "END", "Proportion_of_Unique_Haplotypes"])

    top, num_windows_ls = _avg_from_windows_within_coordinates(data, breakpoint_coordinates, "Proportion_of_Unique_Haplotypes")
    
    null = _permute_data(data, "Proportion_of_Unique_Haplotypes", num_windows_ls, num_permutations)

    z, p = _ztest(top, null)

    _plot_z_test(null, top, "Proportion_of_Unique_Haplotypes", plottitle, z, p)





#ground truth inversion breakpoints
#all breakpoints ###Note: switch to 10Kb on either side of the breakpoint
breakpoint_coords = breakpoint_coordinates_from_bed("inv_SD_bed_files/porubsky_inversions.bed", 10000)

#breakpoints for inversions > 1Mb
large_breakpoint_coords = breakpoint_coordinates_from_bed("inv_SD_bed_files/porubsky_large_inversions.bed", 10000)

#ground truth SD coords ###Note: we are only taking the merged regions and not extending out
SD_coords = SD_coordinates_from_bed("inv_SD_bed_files/vollger_hg38_SDs.merged.bed", 0)

#ground truth SD coords ###Note: we are only taking the merged regions and not extending out
large_SD_coords = SD_coordinates_from_bed("inv_SD_bed_files/large_vollger_hg38_SDs.merged.bed", 0)



# #permutation on inversion breakpoints vs genome
# run_ztest_on_coordinates(
#     "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph", 
#     breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nInversion Breakpoints vs. Whole Genome\n(10,000 SNP windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph", 
#     breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nInversion Breakpoints vs. Whole Genome\n(100Kbp windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph", 
#     breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nInversion Breakpoints vs. Whole Genome\n(1000 SNP windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph", 
#     breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nInversion Breakpoints vs. Whole Genome\n(10Kbp windows)")





# #permutation on breakpoints of large inversions vs genome
# run_ztest_on_coordinates(
#     "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph", 
#     large_breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nLarge Inversion (>1Mb) Breakpoints vs. Whole Genome\n(10,000 SNP windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph", 
#     large_breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nLarge Inversion (>1Mb) Breakpoints vs. Whole Genome\n(100Kbp windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph", 
#     large_breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nLarge Inversion (>1Mb) Breakpoints vs. Whole Genome\n(1000 SNP windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph", 
#     large_breakpoint_coords, 10000, 
#     "Unique Haplotype Proportion\nLarge Inversion (>1Mb) Breakpoints vs. Whole Genome\n(10Kbp windows)")



# #permutation of SDs vs genome
# run_ztest_on_coordinates(
#     "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph", 
#     SD_coords, 10000, 
#     "Unique Haplotype Proportion\nSegmental Duplication Regions vs. Whole Genome\n(10,000 SNP windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph", 
#     SD_coords, 10000, 
#     "Unique Haplotype Proportion\nSegmental Duplication Regions vs. Whole Genome\n(100Kbp windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph", 
#     SD_coords, 10000, 
#     "Unique Haplotype Proportion\nSegmental Duplication Regions vs. Whole Genome\n(1000 SNP windows)")

# run_ztest_on_coordinates(
#     "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph", 
#     SD_coords, 10000, 
#     "Unique Haplotype Proportion\nSegmental Duplication Regions vs. Whole Genome\n(10Kbp windows)")








#permutation of SDs vs genome
run_ztest_on_coordinates(
    "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph", 
    large_SD_coords, 10000, 
    "Unique Haplotype Proportion\nLarge (>100Kb) Segmental Dup Regions vs. Whole Genome\n(10,000 SNP windows)")

run_ztest_on_coordinates(
    "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph", 
    large_SD_coords, 10000, 
    "Unique Haplotype Proportion\nLarge (>100Kb) Segmental Dup Regions vs. Whole Genome\n(100Kbp windows)")

run_ztest_on_coordinates(
    "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph", 
    large_SD_coords, 10000, 
    "Unique Haplotype Proportion\nLarge (>100Kb) Segmental Dup Regions vs. Whole Genome\n(1000 SNP windows)")

run_ztest_on_coordinates(
    "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph", 
    large_SD_coords, 10000, 
    "Unique Haplotype Proportion\nLarge (>100Kb) Segmental Dup Regions vs. Whole Genome\n(10Kbp windows)")
