import os
import pandas as pd
import numpy as np
import pybedtools
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt


###n HELPER FUNCTIONS ###

def _load_bed_by_chrom(bedfile):
    """
    Helper function loads a bedfile of spans as a dictionary of {chromosome:list_of_spans}.
    Where each list contains tuples of (CHROM, START, END). Each list is sorted by start coordinate.
    e.g. {"chr1":[("chr1", 1000, 2000), ("chr1", 5000, 6000)], "chr3":[("chr3", 200, 400), ("chr3", 501, 601)]}

    :param bedfile: filename for bedfile of spans
    :type bedfile: str

    :returns: dictionary that lists spans for each chromosome
    :rtype: dict
    """

    #make empty dictionary to fill
    chrom_to_spans_dict = {}

    #open bedfile and load dictionary with peaks
    with open(bedfile, "r") as peaks_file:
        for line in peaks_file:
            line_list = line.strip().split("\t")
            chrom = line_list[0]
            start = int(line_list[1])
            end = int(line_list[2])
            peak = (chrom, start, end)

            if chrom in chrom_to_spans_dict:
                chrom_to_spans_dict[chrom].append(peak)
            else:
                chrom_to_spans_dict[chrom] = [peak]

    #sort lists
    for l in chrom_to_spans_dict.values():
        l.sort()

    return chrom_to_spans_dict


def _map_peaks_to_invbrks(peaks_dict, invs_dict):
    """
    Helper function that maps peak spans to inversions if they straddle an inversion's breakpoint

    :param peaks_dict: dictionary from `_load_bed_by_chrom` that lists peak spans for each chromosome
    :param invs_dict: dictionary from `_load_bed_by_chrom` that lists inversion spans for each chromosome

    :type peaks_dict: dict
    :type invs_dict: dict

    :returns: dictionary of {peak_tuple:list_of_inversion_tuples}, e.g. {('chr1', 100, 2000): [('chr1', 500, 6000), ('chr1', 70, 105)]}
    :rtype: dict
    """

    #make empty dictionary to fill
    peaks_to_invs_dict = {}

    #get list of chromosome names to run mapping chromosome by chromsome
    chroms = list(set(peaks_dict.keys()))

    #iterate through chromosomes
    for chrom in chroms:
        #grab lists of peaks and inversions
        chrom_peaks = peaks_dict[chrom]
        chrom_invs = invs_dict[chrom]
        
        #iterate through peaks
        for peak in chrom_peaks:
            #add key to dict with empty list
            peaks_to_invs_dict[peak] = []
            #grab peak coordinates
            peak_start = peak[1]
            peak_end = peak[2]

            #iterate through inversions and find mapping
            for inv in chrom_invs:
                #grab inversion coordinates
                inv_start = inv[1]
                inv_end = inv[2]

                #find if peak straddles inversion breakpoint
                if ((inv_start >= peak_start) and (inv_start <= peak_end)) or ((inv_end >= peak_start) and (inv_end <= peak_end)):
                    peaks_to_invs_dict[peak].append(inv)

    return peaks_to_invs_dict


def _list_of_peak_pairs(chrom_to_peaks_dict):
    """
    Helper function takes in the dictionary of peak span lists mapped to chromosomes
    and returns a list of peak span pairs for single peaks, adjacent peaks, and peaks separated by 1 peak

    e.g. {"chr1": [('chr1', 10, 20), ('chr1', 30, 40), ('chr1', 50, 60), ('chr1', 70, 80)], "chr2": [('chr2', 90, 100), ('chr2', 110, 120)]}

    becomes: 
    [(('chr1', 10, 20),('chr1', 10, 20)), (('chr1', 30, 40),('chr1', 30, 40)), ... ,
    (('chr1', 10, 20),('chr1', 30, 40)), (('chr1', 30, 40),('chr1', 50, 60)), (('chr1', 50, 60),('chr1', 70, 80)),
    (('chr1', 10, 20),('chr1', 50, 60)), (('chr1', 30, 40),('chr1', 70, 80)), (('chr2', 90, 100),('chr2', 110, 120))]

    :param chrom_to_peaks_dict: dictionary of {chromosome:list_pf_peak_spans} from `_load_bed_by_chrom`
    :type chrom_to_peaks_dict: dictionary

    :returns: list of tuples of adjacent and separated-by-1 peak pairs
    :rtype: list
    """

    #get list of chromosome names to grab pairs chromosome by chromsome
    chroms = list(set(chrom_to_peaks_dict.keys()))

    #empty list to fill
    peak_pairs = []

    #iterate through chromsomes
    for chrom in chroms:
        #grab list of peaks on chromsome
        list_of_peaks = chrom_to_peaks_dict[chrom]
        #grab list of peak pairs for chromosome
        pairs = list(zip(list_of_peaks, list_of_peaks)) + list(zip(list_of_peaks, list_of_peaks[1:])) + list(zip(list_of_peaks, list_of_peaks[2:]))
        #append to peak pairs list
        peak_pairs += pairs

    return peak_pairs


def _peak_pair_surround_inversion(peak_to_invs_dict, peak_pairs_list):
    """
    Helper function that takes in the dictionary of peak spans mapped to inversion spans whose
    inversion breakpoints they straddle and a list of peak pairs and determines if a peak of peaks
    surrounds an inversion e.g. peak1 and peak2 straddle both breakpoints of an inversion. Note that
    peaks do not overlap, thus an inversion breakpoint cannot be straddled by multiple peaks. So,
    whether a peak pair surrounds an inversion can be determined by whether they share any element in 
    their associated lists of inversion spans from the peak to inverion dict.

    :param peak_to_invs_dict: dictionary of {peak_tuple:list_of_inversion_tuples} from `_map_peaks_to_invbrks`
    :param peak_pairs_list: list of tuples of adjacent and separated by 1 peak pairs from `_list_of_peak_pairs`

    :type peak_to_invs_dict: dict
    :type peak_pairs_list: list

    :return: dataframe with columns made from 1) the list of peak pairs and 2) a list of boolean values (1 or 0) of whether or not the pair surrounds an inversion
    :rtype: pandas.DataFrame
    """

    bool_list = []
    flag_list = []

    #iterate through peak pairs
    for peak1, peak2 in peak_pairs_list:

        if peak1 == peak2:
            flag_list.append("SELF PAIR")
        else:
            flag_list.append("NONSELF PAIR")

        #grab associated lists of inversion
        peak1_invs = peak_to_invs_dict[peak1]
        peak2_invs = peak_to_invs_dict[peak2]

        #does peak pair surround inversion
        bool_val = 0

        #if a peak does not straddle any inversions, then the pair automatically does not surround an inversion
        if len(peak1_invs) == 0 or len(peak2_invs) == 0:
            bool_val = 0
        else:
            #iterate through inversion list1 and check for any identical spans in inversion list2
            for inv in peak1_invs:
                if inv in peak2_invs:
                    bool_val = 1
                    break
                else:
                    bool_val = 0
        
        #append answer to list
        bool_list.append(bool_val)

    #construct df
    df = pd.DataFrame({"peak_pair":peak_pairs_list, "surround_inversion":bool_list, "peak_flag":flag_list})

    return df

# # # # # # # # # # # # #

def _delete_BedTool(BedTool_obj):
    """
    Helper function for deleting the temperorary file for a BedTool object

    :param BedTool_obj: BedTool object instance
    :type BedTool_obj: pybedtools.bedtool.BedTool
    """

    fname = BedTool_obj.fn

    if os.path.exists(fname):
        os.remove(fname)


def _filter_SD_bed(bedfile):
    """
    Helper function takes in the SD bedfile from Vollger et al. 2022 (https://zenodo.org/records/5502036),
    loads it in as a dataframe, filters for SD pairs on the same chromosome and opposite strands, 
    subsets necessary columns, and returns the reindexed dataframe

    :param bedfile: hg38.chr_only.SDs.bed file path
    :type bedfile: str

    :returns: dataframe of filtered bedfile
    :rtype: pandas.DataFrame
    """

    #opening data
    df = pd.read_csv(bedfile, sep="\t")

    #renaming column 1 to get rid of "#"
    df.rename(columns={"#chr1": "chr1"}, inplace=True)

    #filter for SD pairs on same chrom
    df = df[df["chr1"] == df["chr2"]]

    #filter for SD pairs on opposite strand
    df = df[df["strand1"] != df["strand2"]]

    #subset necessary columns only
    df = df[["chr1", "start1", "end1", "strand1", "chr2", "start2", "end2", "strand2", "matchB"]]

    #reset indeces
    df.reset_index(drop=True, inplace=True)
    
    return df


def _SDdf_to_BedTools(SDdf):
    """
    Helper function that converts the sets of SDs from the SD dataframe returned by `_filter_SD_bed`
    and converts them to a BedTool objects so that they can be used by the pybedtools methods. This
    function returns a BedTool of SDs for set 1 (columns: "chr1", "start1", "end1", strand1")
    and for set 2 (columns: "chr2", "start2", "end2", strand2")

    :param SDdf: dataframe from `_filter_SD_bed` that loads in hg38.chr_only.SDs.bed from Vollger 2022
    :type SDdf: pandas.DataFrame

    :returns: BedTool object of set1 SDs and BedTool object of set2 SDs 
    :rtype: pybedtools.bedtool.BedTool, pybedtools.bedtool.BedTool
    """

    #splitting df into set1 and set2
    df1 = SDdf[["chr1", "start1", "end1", "strand1"]]
    df2 = SDdf[["chr2", "start2", "end2", "strand2"]]

    #converting dfs to strings
    SD1_str = df1.to_string(index=False, header=False)
    SD2_str = df2.to_string(index=False, header=False)

    #converting strings to BedTools
    SD1_bed = pybedtools.BedTool(SD1_str, from_string=True)
    SD2_bed = pybedtools.BedTool(SD2_str, from_string=True)

    return SD1_bed, SD2_bed


def _bed_spans_from_peak_pair(peakpair):
    """
    Helper function that takes a peak pair and returns a BedTool object for each peak.

    :param peakpair: pair of peaks stored in a tuple, e.g. (('chr1', 2000, 3000), ('chr1', 5000, 6000))
    :type peakpair: tuple

    :returns: BedTool object of peak 1 and BedTool object of peak 2
    :rtype: pybedtools.bedtool.BedTool, pybedtools.bedtool.BedTool
    """

    #separate peak pair
    peak1_tuple = peakpair[0]
    peak2_tuple = peakpair[1]

    #generate strings for each peak
    peak1_str = f"{peak1_tuple[0]}\t{peak1_tuple[1]}\t{peak1_tuple[2]}"
    peak2_str = f"{peak2_tuple[0]}\t{peak2_tuple[1]}\t{peak2_tuple[2]}"

    #converting to BedTool objects
    peak1_bed = pybedtools.BedTool(peak1_str, from_string=True)
    peak2_bed = pybedtools.BedTool(peak2_str, from_string=True)

    return peak1_bed, peak2_bed


def _find_bedintersect_idxs(intersectbed):
    """
    Helper function loops through a BedTool made from `.intersect` with `loj=True` and returns a
    list of row indeces for row where there was an intersection

    :param intersectbed: BedTool of intersects made with left outer join
    :type intersectbed: pybedtools.bedtool.BedTool

    :returns: list of index values
    :rtype: list
    """

    #empty list to fill with idx values
    idxs = []

    #iterate through intersects and save idx values where intersects occur
    for idx, intersect_interval in enumerate(intersectbed):
        if "-1" not in intersect_interval.fields:
            idxs.append(idx)

    return idxs


def _find_idxs_in_both_sets(idxs1, idxs2):
    """
    Helper function that takes 2 sets of index values adn returns a list of index values
    present in both sets (intersection)

    :param idxs1: list of index values
    :param idxs2: list of index values

    :type idxs1: list
    :type idxs2: list

    :returns: list of index values present in both sets
    :rtype: list
    """

    #converting lists to sets
    set1 = set(idxs1)
    set2 = set(idxs2)

    #finding intersection
    idxs1_2 = list(set1.intersection(set2))

    return idxs1_2


def _find_peak_intersect_indeces(SD1bed, SD2bed, peak1bed, peak2bed):
    """
    Helper function uses pybedtools to find the intersections of SD1 and peak1 and intersections of 
    SD2 and peak2 and note the index positions in the bed dataframe of these intersections. Then the
    function returns a list of indeces present in both intersection sets which give the index positions
    of SDs that map to peak1 and have inverted homologous mates that map to the peak2.

    :param SD1bed: BedTool of SDs for set 1 made with `_SDdf_to_BedTools`
    :param SD2bed: BedTool of SDs for set 2 made with `_SDdf_to_BedTools`
    :param brk1bed: BedTool object of peak 1 made with `_bed_spans_from_peak_pair`
    :param brk2bed: BedTool object of peak 2 made with `_bed_spans_from_peak_pair`

    :type SD1bed: pybedtools.bedtool.BedTool
    :type SD2bed: pybedtools.bedtool.BedTool
    :type brk1bed: pybedtools.bedtool.BedTool
    :type brk2bed: pybedtools.bedtool.BedTool

    :returns: list of index values of rows in the SD dataframe made with `_filter_SD_bed` in which an SD in its inverted homolog mate map to a peak pair
    :rtype: list
    """

    #finding SD intersections with breakpoint 1
    SD_to_peak1 = SD1bed.intersect(peak1bed, loj=True)

    #finding SD intersections with breakpoint 2
    SD_to_peak2 = SD2bed.intersect(peak2bed, loj=True)

    #find idxs for both sets
    sd1_idxs = _find_bedintersect_idxs(SD_to_peak1)
    sd2_idxs = _find_bedintersect_idxs(SD_to_peak2)

    #find idx values present in both sets (intersection)
    sd_intersects_idxs = _find_idxs_in_both_sets(sd1_idxs, sd2_idxs)

    #delete temporary bedfiles
    _delete_BedTool(SD_to_peak1)
    _delete_BedTool(SD_to_peak2)

    return sd_intersects_idxs


def _count_span_results(sddf, intersectidxs):
    """
    Helper function that takes in the filtered SD dataframe and the list of indeces for
    rows of the dataframe where the peak spans map to inverted SD mates and 
    returns the total number of inverted homologous base pairs

    :param sddf: dataframe of filtered bedfile from `_filter_SD_bed`
    :param intersectidxs: list of index values of rows in the SD dataframe from `_find_peak_intersect_indeces`

    :type sddf: pandas.DataFrame
    :type intersectidxs: list

    :returns: number of inverted homologous bps
    :rtype: int
    """

    #counting number of inverted homologus bps
    if len(intersectidxs) == 0:
        num_inv_bps = 0
    else:
        filtered_df = sddf.iloc[intersectidxs]
        num_inv_bps = int(filtered_df["matchB"].sum())

    return num_inv_bps




### HIGHER LEVEL FUNCTIONS ###

def generate_peak_pair_to_inv_df(peaks_bed, inversions_bed):
    """
    Function that takes in a bedfile of peaks and a bed file of inversions and returns a 
    dataframe that maps peak pairs to a boolean value (0 or 1) of whether the pair surrounds an
    inversions. See `_load_bed_by_chrom`, `_map_peaks_to_invbrks`, `_list_of_peak_pairs`,
    and `_peak_pair_surround_inversion` for details.

    :param peaks_bed: filename for bedfile of peak spans
    :param inversions_bed: filename for bedfile of inversion spans

    :type peaks_bed: str
    :type inversions_bed: str

    :returns: dataframe with columns made from 1) the list of peak pairs and 2) a list of boolean values (1 or 0) of whether or not the pair surrounds an inversion
    :rtype: pandas.DataFrame
    """

    #generate by-chromosome dictionaries for peaks and inversions
    chrom_to_peaks = _load_bed_by_chrom(peaks_bed)
    chrom_to_invs = _load_bed_by_chrom(inversions_bed)

    #map peaks to inversions they straddle the breakpoints of
    mapping_peaks_to_invs = _map_peaks_to_invbrks(chrom_to_peaks, chrom_to_invs)

    #generating list of peak pairs
    peak_pairs = _list_of_peak_pairs(chrom_to_peaks)

    #find peaks pairs that surround inversion
    peak_pairs_to_invs = _peak_pair_surround_inversion(mapping_peaks_to_invs, peak_pairs)

    return peak_pairs_to_invs



def count_inverted_bps_for_peak_pairs(peak_pairs_df, SDbed):
    """
    Function that takes in the dataframe from `generate_peak_pair_to_inv_df` and for each 
    peak pair finds the total number of inverted homologous base pairs shared across the 
    2 peaks using the SD bedfile from Vollger et al. 2022 (https://zenodo.org/records/5502036).
    These values added to the dataframe as an additional column.

    :param peak_pairs_df: dataframe with column of peak pairs and column of boolean values (1 or 0) of whether or not the pair surrounds an inversion
    :param SDbed: hg38.chr_only.SDs.bed file path

    :type peak_pairs_df: pandas.DataFrame
    :type SDbed: str

    :returns: same as input dataframe but with an additional column with the number of inverted homologous bps between peak1 and peak2
    :rtype: pandas.DataFrame
    """

    #filter SD bedfile for SDs on the same chromosome and opposite strands
    sd_df = _filter_SD_bed(SDbed)

    #create BedTool objects for each set of SDs
    sd1_bed, sd2_bed = _SDdf_to_BedTools(sd_df)

    #grab list of peak pairs
    peak_pairs_list = peak_pairs_df["peak_pair"].to_list()

    #empty list to fill
    invdup_bps = []

    #iterate through peak pairs and count inverted homology
    for pair in peak_pairs_list:

        #convert peaks to BedTool objects
        peak1_bed, peak2_bed = _bed_spans_from_peak_pair(pair)
        #find inverted SDs that map to both peaks
        mapping_SD_idxs = _find_peak_intersect_indeces(sd1_bed, sd2_bed, peak1_bed, peak2_bed)
        #count number of inverted homooglous bps between the two peaks
        num_bp = _count_span_results(sd_df, mapping_SD_idxs)
        invdup_bps.append(num_bp)

        #delete temporary peak bedfiles
        _delete_BedTool(peak1_bed)
        _delete_BedTool(peak2_bed)

    #delete temporary SD bedfiles
    _delete_BedTool(sd1_bed)
    _delete_BedTool(sd2_bed)

    #add new column to dataframe
    peak_pairs_df["num_invhom_bps"] = invdup_bps

    return peak_pairs_df



def ROC_invhombp_to_invs(peakspairs_inv_bp_df, label):
    """
    Function plots the Receiver Operating Characteristic of the number of inverted homologous bps
    to predict whether a pair of peaks surrounds an inversion.

    :param peakspairs_inv_bp_df: dataframe from `count_inverted_bps_for_peak_pairs`
    :param label: string added to beginning of filename, this is handled by the main functions

    :type peakspairs_inv_bp_df: pandas.DataFrame
    :type label: str
    """

    #classifier parameter
    vals = peakspairs_inv_bp_df["num_invhom_bps"].values

    #groundtruth classification
    truth = peakspairs_inv_bp_df["surround_inversion"].values

    #run ROC
    fpr, tpr, thresholds = roc_curve(truth, vals)
    roc_auc = auc(fpr, tpr)

    #finding Thresholds
    thresh_idx = np.argmin(np.abs(fpr - 0.05))
    fpr_5 = fpr[thresh_idx]
    tpr_5 = tpr[thresh_idx]
    thresh_5 = thresholds[thresh_idx]

    optimal_idx = np.argmax(tpr - fpr)
    fpr_op = fpr[optimal_idx]
    tpr_op = tpr[optimal_idx]
    thresh_op = thresholds[optimal_idx]

    print(label)
    print(f"Threshold with {round(fpr_5, 4)} False Positive Rate and {round(tpr_5, 4)} True Positive Rate is: {int(thresh_5)} bps of inverted homology")
    print(f"Optimal Threshold with {round(fpr_op, 4)} False Positive Rate and {round(tpr_op, 4)} True Positive Rate is: {int(thresh_op)} bps of inverted homology")


    #plotting
    title = f"ROC of # Inverted Homolgous BPs on Identifying\nHigh Hapcount Peaks that Span Inversions"
    plt.figure(figsize=(6,6))
    plt.plot(fpr, tpr, label=f"AUC = {round(roc_auc, 3)}")
    plt.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.7, label="No Skill")
    plt.axvline(x=0.05, color="black", linewidth=0.7, label="5% FPR cutoff")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    plt.savefig(f"{label}{title}_ROC_curve.png")
    plt.clf()



def FDR_invhombp_to_invs(peakspairs_inv_bp_df, label):
    """
    Function plots the TPR vs. FDR and FFR vs. Threshold of the number of 
    inverted homologous bps to predict whether a pair of peaks surrounds an inversion.

    :param peakspairs_inv_bp_df: dataframe from `count_inverted_bps_for_peak_pairs`
    :param label: string added to beginning of filename, this is handled by the main functions

    :type peakspairs_inv_bp_df: pandas.DataFrame
    :type label: str
    """

    #get list of thresholds
    threshs = list(set(peakspairs_inv_bp_df["num_invhom_bps"].to_list()))
    threshs.sort(reverse=True)

    #find number of true positives
    total_tp = peakspairs_inv_bp_df["surround_inversion"].sum()
    total_tn = len(peakspairs_inv_bp_df) - peakspairs_inv_bp_df["surround_inversion"].sum()

    #iterate through thresholds
    tpr = []
    fpr = []
    fdr = []
    for thresh in threshs:
        filtered_df = peakspairs_inv_bp_df[peakspairs_inv_bp_df["num_invhom_bps"] >= thresh]
        
        tp = filtered_df["surround_inversion"].sum()
        fp = len(filtered_df) - filtered_df["surround_inversion"].sum()
        fp_tp = len(filtered_df)

        tpr.append(tp / total_tp)
        fpr.append(fp / total_tn)
        fdr.append(fp / fp_tp)

    #plot
    # plt.figure(figsize=(6,6))
    # plt.plot(fpr, tpr)
    # plt.xlabel("False Positive Rate")
    # plt.ylabel("True Positive Rate")
    # plt.show()

    plt.figure(figsize=(6,6))
    plt.plot(tpr, fdr)
    plt.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.7, label="FDR = TPR")
    plt.ylabel("False Discovery Rate")
    plt.xlabel("True Positive Rate")
    plt.title("False Discovery Rate vs. True Positive Rate\nof # Inverted Homolgous BPs on Identifying Inversions")
    plt.tight_layout()
    
    plt.savefig(f"{label}False Discovery Rate vs. True Positive Rate of # Inverted Homolgous BPs on Identifying Inversions_ROC_curve.png")
    plt.clf()

    plt.figure(figsize=(8,6))
    plt.plot(threshs, fdr)
    plt.ticklabel_format(style='plain', axis='x')
    plt.xlabel("Threshold (bp of inverted homology)")
    plt.ylabel("False Discovery Rate")
    plt.title("False Discovery Rate vs. # Inverted Homolgous BPs as\na Threshold for Identifying Inversions")
    plt.tight_layout()
    
    plt.savefig(f"{label}False Discovery Rate vs. # Inverted Homolgous BPs as\na Threshold for Identifying Inversions_ROC_curve.png")
    plt.clf()












### MAIN FUNCTION ###
def run_program(peakbed, invbed, sdbed):
    """
    Main Function that runs the program to plot a ROC curve of the number of inverted
    homologous bps between a pair of peak to predict an invertion between those peaks.

    :param peakbed: filename for bedfile of peak spans
    :param invbed: filename for bedfile of inversion spans
    :param sdbed: hg38.chr_only.SDs.bed file path

    :type peakbed: str
    :type invbed: str
    :type sdbed: str
    """

    #generate dataframe of peak pairs and whether they surround an inversion
    peaks_inv_df = generate_peak_pair_to_inv_df(peakbed, invbed)

    #add column of count of inverted homoglous bps between peaks
    pairs_inv_bp_df = count_inverted_bps_for_peak_pairs(peaks_inv_df, sdbed)

    #save csv
    pairs_inv_bp_df.to_csv("peakpairs_surroundInvs_invhombp.csv", index=False)

    #filtering for self pairs only and non-self pairs only
    selfpairs_df = pairs_inv_bp_df[pairs_inv_bp_df["peak_flag"] == "SELF PAIR"]
    nonselfpairs_df = pairs_inv_bp_df[pairs_inv_bp_df["peak_flag"] == "NONSELF PAIR"]

    #plot ROC
    ROC_invhombp_to_invs(pairs_inv_bp_df, "")
    #plot FDR
    FDR_invhombp_to_invs(pairs_inv_bp_df, "")

    #plot ROC
    ROC_invhombp_to_invs(selfpairs_df, "Self Pairs Only_")
    #plot FDR
    FDR_invhombp_to_invs(selfpairs_df, "Self Pairs Only_")

    #plot ROC
    ROC_invhombp_to_invs(nonselfpairs_df, "NonSelf Pairs Only_")
    #plot FDR
    FDR_invhombp_to_invs(nonselfpairs_df, "NonSelf Pairs Only_")

    #cleanup any temp files leftover
    pybedtools.cleanup(remove_all=True)




### RUN THE PROGRAM ###

def main():
    run_program("lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed", "porubsky_inversions.bed", "hg38.chr_only.SDs.bed")
    run_program("SD_density_peaks.bed", "porubsky_inversions.bed", "hg38.chr_only.SDs.bed")

if __name__ == "__main__":
    main()
    #df = pd.read_csv("peakpairs_surroundInvs_invhombp.csv")
    #FDR_invhombp_to_invs(df)
