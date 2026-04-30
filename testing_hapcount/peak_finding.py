import pandas as pd
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.signal import find_peaks


def _chrom_sorter(chrom_name):
    """
    Helper function used as key argument in `.sort()` to sort chromosome names by number.
    """

    chrom_num = int(chrom_name.replace("chr", ""))

    return chrom_num


def _get_chromosomes_in_order(chrom_column):
    """
    Helper function used ot get a list of chromosome names in order to be used for filtering later.

    :param chrom_column: CHROM column of dataframe 
    :type chrom_column: pandas.core.series.Series

    :returns: list of chromosome names in order
    :rtype: list
    """

    #find chromosomes
    temp_chroms = set(chrom_column)

    #sort chromosomes
    autosomes = [c for c in temp_chroms if c not in ("chrX", "chrY")]
    autosomes.sort(key=_chrom_sorter)
    chroms = autosomes
    if "chrX" in temp_chroms:
        chroms.append("chrX")
    if "chrY" in temp_chroms:
        chroms.append("chrY")

    return chroms



def _identify_spans_from_lowess(boolean_above_threshold):
    """
    Helper function that returns the index positions for the starts and stops for spans of data above the threshold.
    This helper is used in the `lowess_thresholding` function. It returns a list of tuples with start and stop indeces.
    These indeces can then be mapped back to the chrom_bedgraph_df to find the start and stop genomic coordinates.

    :param boolean_above_threshold: array of 0s and 1s representing whether data in that index of the chrom_bedgraph_df is above the threshold
    :type boolean_above_threshold: numpy.ndarray

    :returns: list of tuples of 2 integers with a tuple for span above the threshold and the integers being the start and stop indeces of the chrom_bedgraph_df for the span
    :rtype: [(int, int)]
    """

    #finding runs above the threshold, padding ends with 0
    diff = np.diff(boolean_above_threshold.astype(int), prepend=0, append=0)

    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    starts_and_stops = list(zip(starts, ends))

    return starts_and_stops


def lowess_thresholding(chrom_bedgraph_df, percentile_threshold, fraction_param):
    """
    Function that identifies peaks of windowed genomic data by first lowess smoothing the data and then identifying 
    spans above a specified percentile of the smoothed data. Note that this function should only be applied one chromosome at a time.

    :param chrom_bedgraph_df: dataframe made from the input bedgraph file with the column names "CHROM", "START", "END", and "DATA"
    :param percentile_threshold: percentlile (given as a percent value) of the smoothed data above which to idnetify peaks
    :param fraction_param: frac parameter to be used in lowess smoothing, 0.005 is recommended. See `statsmodels.nonparametric.smoothers_lowess.lowess` docs

    :type chrom_bedgraph_df: DataFrame
    :type percentile_threshold: int or float
    :type fraction_param: float
    :type prominence_param: float

    :returns: list of tuples representing the genomic spans for each identified peak
    :rtype: [(str, int, int)]
    """

    print("- finding peaks using lowess smoothed data")

    #finding chromosome
    chrom = chrom_bedgraph_df["CHROM"][0]

    #lowess smoothing
    lowess_smooth = lowess(chrom_bedgraph_df["DATA"], chrom_bedgraph_df["START"], frac=fraction_param, return_sorted=False)

    #finding percentile threshold
    pct_thresh = np.percentile(lowess_smooth, percentile_threshold)

    #finding bool list for lowess data >= threshold
    bool_above_thresh = lowess_smooth >= pct_thresh

    #finding start and stop idxs
    start_and_stop_idxs = _identify_spans_from_lowess(bool_above_thresh)


    #constructing list of genomic spans
    spans = []
    for start_idx, stop_idx in start_and_stop_idxs:
        c = chrom
        s = chrom_bedgraph_df.at[start_idx, "START"]
        e = chrom_bedgraph_df.at[stop_idx-1, "END"] #-1 to avoid off by one error at end of chromosome and to not include window outside of span

        spans.append((c, s, e))

    return spans




def height_and_prominence_thresholding(chrom_bedgraph_df, percentile_threshold, prominence_param):
    """
    Function that identifies peaks of windowed genomic data by using the `scipy.signal.find_peaks` algorithm.
    A percentile is used to idenify peaks above this threshold and of these peaks are kept if their prominence (see scipy docs)
    is above a given value. Note that this function should only be applied one chromosome at a time.

    :param chrom_bedgraph_df: dataframe made from the input bedgraph file with the column names "CHROM", "START", "END", and "DATA"
    :param percentile_threshold: percentlile of the (not smoothed) data above which to idnetify peaks
    :param prominence_param: promience parameter to be using in peak finding, 0.25 is recommended. See `scipy.signal.find_peaks` docs

    :type chrom_bedgraph_df: DataFrame
    :type percentile_threshold: int or float

    :returns: list of tuples representing the genomic spans for each identified peak
    :rtype: [(str, int, int)]
    """

    print("- finding peaks using find_peaks()")

    #finding chromosome
    chrom = chrom_bedgraph_df["CHROM"][0]

    #finding percentile threshold
    pct_thresh = np.percentile(chrom_bedgraph_df["DATA"], percentile_threshold)

    #finding peaks
    peaks, properties = find_peaks(chrom_bedgraph_df["DATA"], height=pct_thresh, prominence=prominence_param, width=0)

    #constructing list of genomic spans
    spans = []
    for idx in range(len(properties["left_ips"])):

        #find row indeces
        left_idx = int(properties["left_ips"][idx])
        right_idx = int(properties["right_ips"][idx])

        #find coordinates
        c = chrom
        s = chrom_bedgraph_df.at[left_idx, "START"]
        e = chrom_bedgraph_df.at[right_idx, "END"]

        #appending
        spans.append((c, s, e))

    return spans



def lowess_height_and_prominence_thresholding(chrom_bedgraph_df, percentile_threshold, fraction_param, prominence_param):
    """
    Function that identifies peaks of windowed genomic data by first lowess smoothing the data and then uses the `scipy.signal.find_peaks` algorithm
    to idenity peaks in the smoothed data with specified percentile and prominence thresholds. Note that this function should only be applied one chromosome at a time.

    :param chrom_bedgraph_df: dataframe made from the input bedgraph file with the column names "CHROM", "START", "END", and "DATA"
    :param percentile_threshold: percentlile of the smoothed data above which to idnetify peaks
    :param fraction_param: frac parameter to be used in lowess smoothing, 0.005 is recommended. See `statsmodels.nonparametric.smoothers_lowess.lowess` docs
    :param prominence_param: promience parameter to be using in peak finding, 0.05 is recommended. See `scipy.signal.find_peaks` docs

    :type chrom_bedgraph_df: DataFrame
    :type percentile_threshold: int or float
    :type fraction_param: float
    :type prominence_param: float

    :returns: list of tuples representing the genomic spans for each identified peak
    :rtype: [(str, int, int)]
    """

    print("- finding peaks using find_peaks() on lowess smoothed data")

    #finding chromosome
    chrom = chrom_bedgraph_df["CHROM"][0]

    #lowess smoothing
    lowess_smooth = lowess(chrom_bedgraph_df["DATA"], chrom_bedgraph_df["START"], frac=fraction_param, return_sorted=False)

    #finding percentile threshold
    pct_thresh = np.percentile(lowess_smooth, percentile_threshold)

    #finding peaks
    peaks, properties = find_peaks(lowess_smooth, height=pct_thresh, prominence=prominence_param, width=0)

    #constructing list of genomic spans
    spans = []
    for idx in range(len(properties["left_ips"])):

        #find row indeces
        left_idx = int(properties["left_ips"][idx])
        right_idx = int(properties["right_ips"][idx])

        #find coordinates
        c = chrom
        s = chrom_bedgraph_df.at[left_idx, "START"]
        e = chrom_bedgraph_df.at[right_idx, "END"]

        #appending
        spans.append((c, s, e))

    return spans








def find_peaks_from_bedgraph_v1(bedgraph):
    """
    Function takes in a bedgraph file and parses it by chromosome by converting it to a pandas DataFrame with columns
    "CHROM", "START", "END", and "DATA" and filtering by chromsome. Each chromosome df is then passed to the peak identifying functions.
    The function outputs a list for each peak idenitfication method used.

    :param bedgraph: filename for the input bedgraph
    :type bedgraph: str

    :returns: 5 lists for each of the peak identification methods used
    :rtype: list, list, list, list, list
    """

    #loading data
    print("LOADING BEDGRAPH FILE...")
    df = pd.read_csv(bedgraph, sep="\t", header=None, names=["CHROM", "START", "END", "DATA"])

    #grabbing chromosomes
    chromosomes = _get_chromosomes_in_order(df["CHROM"])

    #setting empty peak span lists
    lo90 = []
    lo95 = []
    fp90 = []
    fp95 = []
    lofp90 = []

    #iterating through chromosomes
    print("FILTERING BY CHROMOSOME...")
    for chrm in chromosomes:

        #filtering for chromosome
        chrom_df = df[df["CHROM"] == chrm]
        chrom_df.reset_index(inplace=True, drop=True)

        #finding peaks with all 5 methods
        print(f"RUNNING ON {chrm}")
        lo90spans = lowess_thresholding(chrom_df, 90, 0.005)
        lo95spans = lowess_thresholding(chrom_df, 95, 0.005)
        fp90spans = height_and_prominence_thresholding(chrom_df, 90, 0.25)
        fp95spans = height_and_prominence_thresholding(chrom_df, 95, 0.25)
        lofp90spans = lowess_height_and_prominence_thresholding(chrom_df, 90, 0.005, 0.05)

        #concatenating spans to full genome span lists
        lo90 += lo90spans
        lo95 += lo95spans
        fp90 += fp90spans
        fp95 += fp95spans
        lofp90 += lofp90spans

    return lo90, lo95, fp90, fp95, lofp90





def write_span_list_to_bed(span_list, outbed):
    """
    Function takes a list of spans from `find_peaks_from_bedgraph_v1` and writes to an output bed file.
    
    :param span_list: list of genomic spans for each peak
    :param outbed: filename for output bed file

    :type span_list: [(str, int, int)]
    :type outbed: str
    """

    with open(outbed, "a") as output:
        for c, s, e in span_list:
            output.write(f"{c}\t{s}\t{e}\n")





def main(bedgraph):
    """
    Main function that runs the program on an input bedgraph file
    """

    lowess90, lowess95, findpeaks90, findpeaks95, lowessfindpeaks90 = find_peaks_from_bedgraph_v1(bedgraph)

    write_span_list_to_bed(lowess90, f"RAW_lowess90_{bedgraph}".replace(".bedgraph", ".bed"))
    write_span_list_to_bed(lowess95, f"RAW_lowess95_{bedgraph}".replace(".bedgraph", ".bed"))
    write_span_list_to_bed(findpeaks90, f"RAW_findpeaks90_{bedgraph}".replace(".bedgraph", ".bed"))
    write_span_list_to_bed(findpeaks95, f"RAW_findpeaks95_{bedgraph}".replace(".bedgraph", ".bed"))
    write_span_list_to_bed(lowessfindpeaks90, f"RAW_lowessfindpeaks90_{bedgraph}".replace(".bedgraph", ".bed"))
    


main("Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph")
main("Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph")