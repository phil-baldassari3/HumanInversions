import pandas as pd
import numpy as np


def _mark_windows_within_coordinates(df, coordinates, label):
    """
    Filters a given dataframe by the specified genomic coordinates and outputs a list of labels for each window specifying whether they are within given coordinataes

    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param coordinates: list of tuples with regions to include [("chrom", start, end)...]
    :type coordinates: list of tuples of str, int, and int
    :param label: what to label each window that is within coordinates, else the label will be "none"
    :type label: str

    :returns: list of labels to be added to the dataframe as a new column
    :rtype: list
    """

    #setting mask
    keep_rows = np.zeros(len(df), dtype=bool)

    #bool masking rows for windows entirely or partially overlapping with regions to keep
    for chrom, region_start, region_end in coordinates:
        keep_rows |= ((df["CHROM"] == chrom) & (df["START"] < region_end) & (df["END"] > region_start))

    #making list
    labs = []
    for i in keep_rows:
        if i == 0:
            labs.append("none")
        elif i == 1:
            labs.append(label)
        else:
            print("Something is wrong!")

    
    return labs




def _add_label_col(df, coords, colname, lab):
    """
    Adds a column of labels for each window given a list of coordinates, column name, abnd label word
    
    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param coords: list of tuples with regions to include [("chrom", start, end)...]
    :type coords: list of tuples of str, int, and int
    :param colname: name of column to be added to dataframe
    :type colname: str
    :param lab: label to be used to mark windows within coordinates
    :type lab: str

    :returns: dataframe with new label column
    :rtype: DataFrame
    """

    #getting label list
    lab_list = _mark_windows_within_coordinates(df, coords, lab)

    #adding col
    df[colname] = lab_list

    return df



def add_labels_for_inv_categories(tsv, list_of_coord_lists, list_of_inv_categories):
    """
    Labels windows in tsv file for whether they are in inversions, large inversions, inversion breakpoints, and large inversion breakpoints and saves a new tsv.
    Notes that the label word and column added to the table for each category will be the same
    
    :param tsv: tsv file of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type tsv: str
    :param list_of_coord_lists: list of lists of tuples of coordinates to label windows in
    :type list_of_coord_lists: list
    :param list_of_inv_categories: list of inversion categroies (inversion, large_inversion, etc.) this list MUST be in the same order as list_of_coord_lists
    :type list_of_inv_categories: list
    """

    #opening data
    data = pd.read_csv(tsv, sep="\t")

    #adding columns
    for idx, cat in enumerate(list_of_inv_categories):
        newcolname = cat
        insidelabel = cat
        coords_to_use = list_of_coord_lists[idx]

        data = _add_label_col(data, coords_to_use, newcolname, insidelabel)

    #outputting new file
    data.to_csv(tsv.replace(".tsv", "_labeled.tsv"), sep="\t", index=False)




#setting lists
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



coords_lists = [inv_coords, large_inv_coords, breakpoint_coords, large_breakpoint_coords]
categories = ["inversion", "large_inversion", "inv_breakpoint", "large_inv_breakpoint"]


add_labels_for_inv_categories("hapcount_SNPwindow/cleaned_autosomes_SNPwindow1000_SNPstep500_hapcount_snpden.tsv", coords_lists, categories)
add_labels_for_inv_categories("hapcount_BPwindow/cleaned_autosomes_BPwindow3000_BPstep1500_hapcount_snpden.tsv", coords_lists, categories)