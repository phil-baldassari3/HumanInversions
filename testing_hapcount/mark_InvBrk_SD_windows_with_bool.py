import pandas as pd
import numpy as np




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




def _mark_windows_within_coordinates_with_bool(df, coordinates):
    """
    Filters a given dataframe by the specified genomic coordinates and outputs a list of labels for each window specifying whether they are within given coordinataes with a 1 or 0

    :param df: dataframe of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type df: DataFrame
    :param coordinates: list of tuples with regions to include [("chrom", start, end)...]
    :type coordinates: list of tuples of str, int, and int
    :param label: what to label each window that is within coordinates, else the label will be "none"
    :type label: str

    :returns: list 1s and 0s to be added to the dataframe as a new column
    :rtype: list
    """

    #setting mask
    keep_rows = np.zeros(len(df), dtype=bool)

    #bool masking rows for windows entirely or partially overlapping with regions to keep
    for chrom, region_start, region_end in coordinates:
        keep_rows |= ((df["CHROM"] == chrom) & (df["START"] < region_end) & (df["END"] > region_start))

    #making list
    bools = []
    for i in keep_rows:
        if i:
            bools.append(1)
        else:
            bools.append(0)

    
    return bools




def _add_bool_col(df, coords, colname):
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
    bool_list = _mark_windows_within_coordinates_with_bool(df, coords)

    #adding col
    df[colname] = bool_list

    return df



def add_labels_for_categories(bedgraph, list_of_coord_lists, list_of_inv_categories):
    """
    Labels windows in tsv file for whether they are in inversions, large inversions, inversion breakpoints, and large inversion breakpoints and saves a new tsv.
    Notes that the label word and column added to the table for each category will be the same
    
    :param bedgraph: tsv file of windowed average branch length, repeat density, SNP density, and any other relevant parameters
    :type tsv: str
    :param list_of_coord_lists: list of lists of tuples of coordinates to label windows in
    :type list_of_coord_lists: list
    :param list_of_inv_categories: list of inversion categroies (inversion, large_inversion, etc.) this list MUST be in the same order as list_of_coord_lists
    :type list_of_inv_categories: list
    """

    #opening data
    data = pd.read_csv(bedgraph, sep="\t", header=None, names=["CHROM", "START", "END", "Proportion_of_Unique_Haplotypes"])

    #adding columns
    for idx, cat in enumerate(list_of_inv_categories):
        newcolname = cat
        coords_to_use = list_of_coord_lists[idx]

        data = _add_bool_col(data, coords_to_use, newcolname)

    #outputting new file
    data.to_csv(bedgraph.split("/")[-1].replace(".bedgraph", "_boollabeled.tsv"), sep="\t", index=False)







#setting lists
#ground truth sets
#all breakpoints ###Note: switch to 10Kb on either side of the breakpoint
breakpoint_coords = breakpoint_coordinates_from_bed("inv_SD_bed_files/porubsky_inversions.bed", 10000)

#breakpoints for inversions > 1Mb
large_breakpoint_coords = breakpoint_coordinates_from_bed("inv_SD_bed_files/porubsky_large_inversions.bed", 10000)

#ground truth SD coords ###Note: we are only taking the merged regions and not extending out
SD_coords = SD_coordinates_from_bed("inv_SD_bed_files/vollger_hg38_SDs.merged.bed", 0)

#ground truth SD coords ###Note: we are only taking the merged regions and not extending out
large_SD_coords = SD_coordinates_from_bed("inv_SD_bed_files/large_vollger_hg38_SDs.merged.bed", 0)


###Running the program###
coords_lists = [breakpoint_coords, large_breakpoint_coords, SD_coords, large_SD_coords]
categories = ["InvBrk", "LargeInvBrk", "SD", "LargeSD"]


add_labels_for_categories(
    "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph", 
    coords_lists, categories
    )


add_labels_for_categories(
    "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph", 
    coords_lists, categories
    )


add_labels_for_categories(
    "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph", 
    coords_lists, categories
    )


add_labels_for_categories(
    "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph", 
    coords_lists, categories
    )