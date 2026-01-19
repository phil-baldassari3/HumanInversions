import os,sys
import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev
import scipy.stats
import random
import matplotlib.pyplot as plt


random.seed(672)



def _ztest(top_value, null_set):
    """
    Performs a two-tailed Z test between the set of interest and the null set

    :param top_value: integer or float for the average value from the set of interest
    :type top_value: float
    :param null_set: list of integers or floats that of average values in the null set from the permutation
    :type null_set: list of floats
    
    :returns: z-score and p-value
    :rtype: (float, float)
    """

    #to avoid division by zero error but this means that the z test is meaningless
    if stdev(null_set) == 0:
        zscore = 1
        pvalue = 1

    else:   
        zscore = (top_value - mean(null_set)) / stdev(null_set)
        pvalue = 2 * scipy.stats.norm.sf(zscore)

    return zscore, pvalue



def _avg_nth_percentile_of_data(df, pct, pct_param, val_param):
    """
    Filters a given dataframe by the nth percentile of a specified parameter and outputs the average corresponding values from the dataframe
    
    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param pct: percentile to filter dataframe for
    :type pct: float
    :param pct_param: column name for which to filter by percentile on
    :type pct_param: str
    :param val_param: column name for which to grab values after filtering
    :type val_param: str

    :returns: average of values after filtering and the number of windows in the nth percentile
    :rtype: (float, int)
    """

    #finding percentile
    data_to_find_pct = df[pct_param].to_list()
    pct_val = np.percentile(data_to_find_pct, pct)

    #filtering
    filtered_df = df[df[pct_param] >= pct_val]

    #grabbing values
    val_list = filtered_df[val_param].to_list()
    sample_size = len(val_list)

    #averaging values
    avg_val = mean(val_list)

    return avg_val, sample_size


def _avg_from_windows_within_coordinates(df, coordinates, val_param):
    """
    Filters a given dataframe by the specified genomic coordinates and outputs the average of specified values from the dataframe

    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param coordinates: list of tuples with regions to include [("chrom", start, end)...]
    :type coordinates: list of tuples of str, int, and int
    :param val_param: column name for which to grab values after filtering
    :type val_param: str

    :returns: average of values after filtering and the number of windows included
    :rtype: (float, int)
    """

    #setting mask
    keep_rows = np.zeros(len(df), dtype=bool)

    #bool masking rows for windows entirely or partially overlapping with regions to keep
    for chrom, region_start, region_end in coordinates:
        keep_rows |= ((df["CHROM"] == chrom) & (df["START"] < region_end) & (df["END"] > region_start))

    #keeping rows within coordinates
    filtered_df = df[keep_rows]

    #grabbing values
    val_list = filtered_df[val_param].to_list()
    sample_size = len(val_list)

    #averaging values
    avg_val = mean(val_list)

    return avg_val, sample_size



def _permute_data(df, param_to_permute, permutations, n):
    """
    Permutes data a given number of times, each time grabbing a random sample of a specified 
    sample size of windows, finding the average, and appending to a list to generate a null distribution.
    
    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param param_to_permute: column name for values to sample and average
    :param permutations: number of random permutations
    :param n: sample size of windows to sample for each permutation

    :returns: lists of averaged values, one for each permutation thus generating a null distribution
    :rtype: list of floats
    """

    #grabbing data
    vals_to_permute = df[param_to_permute].to_list()

    #empty list for null distribution
    null_dist = []

    for i in range(permutations):
        if i % 100 == 0:
            print(f"Permutation #{i}")
        sampled_vals = random.sample(vals_to_permute, k=n)
        avg_of_sample = mean(sampled_vals)
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
    plt.xlabel(param_name)
    plt.title(plot_title)
  
    plt.axvline(x=val_of_interest, color='r', lw=3)

    plt.text(0.5, 0.8, f"z-score = {round(zscore, 3)}\np-value = {pvalue}", 
             horizontalalignment='center', verticalalignment='bottom', bbox=dict(facecolor='yellow', alpha=0.25),
             transform=plt.gca().transAxes)

    plt.savefig(f"{plot_title}.png", bbox_inches='tight')
    #plt.show()
    plt.clf()





def run_ztest_on_percentile(tsvfile, percentile, param_to_define_percentile, param_to_perform_ztest, num_permutations, plottitle):
    """
    Main function to run permutation Z-test on a specified parameter to test for a significant difference of this parameter
    between its values above the nth percentile of another paramter and its values in the dataset as a whole.
    
    :param tsvfile: file of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type tsvfile: str
    :param percentile: percentile to filter data
    :type percentile: float
    :param param_to_define_percentile: column name for parameter on which to define percentile threshold
    :type param_to_define_percentile: str
    :param param_to_perform_ztest: column name for parameter on which to run the permutation Z-test
    :type param_to_perform_ztest: str
    :param num_permutations: number of permutations
    :type num_permutations: int
    :param plottitle: Title of plot and name of filename
    :type plottitle: str
    """

    data = pd.read_csv(tsvfile, sep="\t")

    top, size = _avg_nth_percentile_of_data(data, percentile, param_to_define_percentile, param_to_perform_ztest)

    null = _permute_data(data, param_to_perform_ztest, num_permutations, size)

    z, p = _ztest(top, null)

    _plot_z_test(null, top, param_to_perform_ztest, plottitle, z, p)




def run_ztest_on_coordinates(tsvfile, coordinates_to_extract, param_to_perform_ztest, num_permutations, plottitle):
    """
    Main function to run permutation Z-test on a specified parameter to test for a significant difference of this parameter
    between its values above the nth percentile of another paramter and its values in the dataset as a whole.
    
    :param tsvfile: file of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type tsvfile: str
    :param coordinates_to_extract: list of tuples with regions to include [("chrom", start, end)...]
    :type coordinates_to_extract: list of tuples of str, int, and int
    :param param_to_perform_ztest: column name for parameter on which to run the permutation Z-test
    :type param_to_perform_ztest: str
    :param num_permutations: number of permutations
    :type num_permutations: int
    :param plottitle: Title of plot and name of filename
    :type plottitle: str
    """

    data = pd.read_csv(tsvfile, sep="\t")

    top, size = _avg_from_windows_within_coordinates(data, coordinates_to_extract, param_to_perform_ztest)

    null = _permute_data(data, param_to_perform_ztest, num_permutations, size)

    z, p = _ztest(top, null)

    _plot_z_test(null, top, param_to_perform_ztest, plottitle, z, p)



# run_ztest_on_percentile("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", 75, "Avg_Branch_Length", "SNP_Density", 1000, "SNP Density of Windows in 75th Branch Length Percentile\nCompared to Whole Genome")
# run_ztest_on_percentile("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", 90, "Avg_Branch_Length", "SNP_Density", 1000, "SNP Density of Windows in 90th Branch Length Percentile\nCompared to Whole Genome")
# run_ztest_on_percentile("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", 98, "Avg_Branch_Length", "SNP_Density", 1000, "SNP Density of Windows in 98th Branch Length Percentile\nCompared to Whole Genome")

# run_ztest_on_percentile("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", 75, "Avg_Branch_Length", "Repeat_Density", 1000, "Repeat Density of Windows in 75th Branch Length Percentile\nCompared to Whole Genome")
# run_ztest_on_percentile("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", 90, "Avg_Branch_Length", "Repeat_Density", 1000, "Repeat Density of Windows in 90th Branch Length Percentile\nCompared to Whole Genome")
# run_ztest_on_percentile("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", 98, "Avg_Branch_Length", "Repeat_Density", 1000, "Repeat Density of Windows in 98th Branch Length Percentile\nCompared to Whole Genome")



#all ground truth inversions
inv_coords = [
    ("chr8", 7064966-50000, 12716088+50000),
    ("chr15", 30077909-50000, 32607507+50000),
    ("chr3", 195622595-50000, 197667517+50000),
    ("chr17", 36357258-50000, 37957978+50000),
    ("chr15", 74060644-50000, 75304798+50000),
    ("chr17", 45495836-50000, 46707123+50000),
    ("chr19", 21647331-50000, 22062459+50000),
    ("chr17", 18597986-50000, 18848096+50000),
    ("chr7", 54223281-50000, 54319128+50000),
    ("chrX", 72996088-50000, 73086935+50000),
    ("chrY", 15845242-50000, 15921740+50000),
    ("chrX", 154338058-50000, 154394093+50000),
    ("chr13", 63716851-50000, 63769457+50000),
    ("chrX", 106266497-50000, 106299864+50000),
    ("chr6", 141866324-50000, 141898714+50000),
    ("chr13", 79819767-50000, 79843107+50000),
    ("chr2", 240676144-50000, 240698757+50000),
    ("chr3", 162807955-50000, 162829867+50000),
    ("chr14", 34540845-50000, 34562270+50000),
    ("chr11", 71564791-50000, 71583752+50000),
    ("chrX", 101597544-50000, 101616245+50000),
    ("chr7", 70955982-50000, 70973901+50000),
    ("chr16", 75205642-50000,75223319+50000),
    ("chr9", 123976373-50000, 123993772+50000),
    ("chr5", 64466001-50000, 64482016+50000),
    ("chr9", 113105263-50000, 113119310+50000),
    ("chr17", 30616735-50000, 30630735+50000),
    ("chr4", 87926012-50000, 87937548+50000)
]

#inversions > 1Mb
large_inv_coords = [
    ("chr8", 7064966-50000, 12716088+50000),
    ("chr15", 30077909-50000, 32607507+50000),
    ("chr3", 195622595-50000, 197667517+50000),
    ("chr17", 36357258-50000, 37957978+50000),
    ("chr15", 74060644-50000, 75304798+50000),
    ("chr17", 45495836-50000, 46707123+50000)
]

#all breakpoints
breakpoint_coords = [
    ("chr8", 7064966-50000, 7064966+50000),
    ("chr15", 30077909-50000, 30077909+50000),
    ("chr3", 195622595-50000, 195622595+50000),
    ("chr17", 36357258-50000, 36357258+50000),
    ("chr15", 74060644-50000, 74060644+50000),
    ("chr17", 45495836-50000, 45495836+50000),
    ("chr19", 21647331-50000, 21647331+50000),
    ("chr17", 18597986-50000, 18597986+50000),
    ("chr7", 54223281-50000, 54223281+50000),
    ("chrX", 72996088-50000, 72996088+50000),
    ("chrY", 15845242-50000, 15845242+50000),
    ("chrX", 154338058-50000, 154338058+50000),
    ("chr13", 63716851-50000, 63716851+50000),
    ("chrX", 106266497-50000, 106266497+50000),
    ("chr6", 141866324-50000, 141866324+50000),
    ("chr13", 79819767-50000, 79819767+50000),
    ("chr2", 240676144-50000, 240676144+50000),
    ("chr3", 162807955-50000, 162807955+50000),
    ("chr14", 34540845-50000, 34540845+50000),
    ("chr11", 71564791-50000, 71564791+50000),
    ("chrX", 101597544-50000, 101597544+50000),
    ("chr7", 70955982-50000, 70955982+50000),
    ("chr16", 75205642-50000,75205642+50000),
    ("chr9", 123976373-50000, 123976373+50000),
    ("chr5", 64466001-50000, 64466001+50000),
    ("chr9", 113105263-50000, 113105263+50000),
    ("chr17", 30616735-50000, 30616735+50000),
    ("chr4", 87926012-50000, 87926012+50000),
    
    ("chr8", 12716088-50000, 12716088+50000),
    ("chr15", 32607507-50000, 32607507+50000),
    ("chr3", 197667517-50000, 197667517+50000),
    ("chr17", 37957978-50000, 37957978+50000),
    ("chr15", 75304798-50000, 75304798+50000),
    ("chr17", 46707123-50000, 46707123+50000),
    ("chr19", 22062459-50000, 22062459+50000),
    ("chr17", 18848096-50000, 18848096+50000),
    ("chr7", 54319128-50000, 54319128+50000),
    ("chrX", 73086935-50000, 73086935+50000),
    ("chrY", 15921740-50000, 15921740+50000),
    ("chrX", 154394093-50000, 154394093+50000),
    ("chr13", 63769457-50000, 63769457+50000),
    ("chrX", 106299864-50000, 106299864+50000),
    ("chr6", 141898714-50000, 141898714+50000),
    ("chr13", 79843107-50000, 79843107+50000),
    ("chr2", 240698757-50000, 240698757+50000),
    ("chr3", 162829867-50000, 162829867+50000),
    ("chr14", 34562270-50000, 34562270+50000),
    ("chr11", 71583752-50000, 71583752+50000),
    ("chrX", 101616245-50000, 101616245+50000),
    ("chr7", 70973901-50000, 70973901+50000),
    ("chr16", 75223319-50000,75223319+50000),
    ("chr9", 123993772-50000, 123993772+50000),
    ("chr5", 64482016-50000, 64482016+50000),
    ("chr9", 113119310-50000, 113119310+50000),
    ("chr17", 30630735-50000, 30630735+50000),
    ("chr4", 87937548-50000, 87937548+50000)
]


#breakpoints for inversions > 1Mb
large_breakpoint_coords = [
    ("chr8", 7064966-50000, 7064966+50000),
    ("chr15", 30077909-50000, 30077909+50000),
    ("chr3", 195622595-50000, 195622595+50000),
    ("chr17", 36357258-50000, 36357258+50000),
    ("chr15", 74060644-50000, 74060644+50000),
    ("chr17", 45495836-50000, 45495836+50000),

    ("chr8", 12716088-50000, 12716088+50000),
    ("chr15", 32607507-50000, 32607507+50000),
    ("chr3", 197667517-50000, 197667517+50000),
    ("chr17", 37957978-50000, 37957978+50000),
    ("chr15", 75304798-50000, 75304798+50000),
    ("chr17", 46707123-50000, 46707123+50000)
]




run_ztest_on_coordinates("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", inv_coords, "Avg_Branch_Length", 1000, "Average Branch Length for Inversion Windows\nCompared to Whole Genome")
run_ztest_on_coordinates("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", large_inv_coords, "Avg_Branch_Length", 1000, "Average Branch Length for Large Inversion (>1Mb) Windows\nCompared to Whole Genome")

run_ztest_on_coordinates("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", breakpoint_coords, "Avg_Branch_Length", 1000, "Average Branch Length for Inversion Breakpoint Windows\nCompared to Whole Genome")
run_ztest_on_coordinates("cleaned_windowed_AvgBranchLen_RepeatDen_SNPden.tsv", large_breakpoint_coords, "Avg_Branch_Length", 1000, "Average Branch Length for Large Inversion (>1Mb) Breakpoint Windows\nCompared to Whole Genome")

