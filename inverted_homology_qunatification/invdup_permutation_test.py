import os
import pandas as pd
import numpy as np
from statistics import mean, stdev
import random
import scipy.stats
import matplotlib.pyplot as plt
from multiprocessing import Pool
import pybedtools
random.seed(999)


### HELPER FUNCTIONS ###

def _generate_seed_value():
    """
    Helper function that returns random values that can be used as seed values 
    for `pybedtools.BedTool.shuffle` calls. Make sure you set a global seed value
    using `random`

    :returns: random number
    :rtype: int
    """
    r = random.randrange(0, 999999999999)
    return r



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



def _load_inversion_bed(invbed):
    """
    Helper function that loads in a bed file of inversion coordinates and returns each inversion
    in a list of tuples where each tuple contains the chrom, start, and end of each span

    :param invbed: bedfile name of inversions
    :type invbed: str

    :returns: list of tuples where each tuple is an inversion and contains (chrom, start, end)
    :rtype: [tuple]
    """

    #list of inversion coordinates to be filled
    inv_list = []

    #opening bed and iterating through lines
    with open(invbed, "r") as bed:
        for line in bed:
            line_list = line.strip().split("\t")
            chrom = line_list[0]
            start = int(line_list[1])
            end = int(line_list[2])

            inv_list.append((chrom, start, end))

    return inv_list



def _brkpnt_bed_spans_from_inv(inv, bp):
    """
    Helper function that takes in a tuple of an inversion coordinate and outputs a BedTool objects
    for each inversion breakpoint coordinate +/- the specified number of base pairs

    e.g. with bp=10 
    chr1    123    456
    becomes
    chr1    113    133 and chr1    446    466

    :param inv: tuple of inversion span: e.g. (chr1, 123, 456)
    :param bp: number of base pairs on either side of each span edge

    :type inv: (str, int, int)
    :type bp: int

    :returns: BedTool object of breakpoint 1 and BedTool object of  breakpoint 2
    :rtype: pybedtools.bedtool.BedTool, pybedtools.bedtool.BedTool
    """

    #getting breakpoint span strings
    chrom = inv[0]
    start1 = max(0, inv[1]-bp)
    end1 = inv[1]+bp
    start2 = max(0, inv[2]-bp)
    end2 = inv[2]+bp

    brk_str1 = f"{chrom}\t{start1}\t{end1}"
    brk_str2 = f"{chrom}\t{start2}\t{end2}"

    #converting to BedTool objects
    brk1_bed = pybedtools.BedTool(brk_str1, from_string=True)
    brk2_bed = pybedtools.BedTool(brk_str2, from_string=True)

    return brk1_bed, brk2_bed



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
         


def _find_brkpnt_intersect_indeces(SD1bed, SD2bed, brk1bed, brk2bed):
    """
    Helper function uses pybedtools to find the intersections of SD1 and brk1 and intersections of 
    SD2 and brk2 and note the index positions in the bed dataframe of these intersections. Then the
    function returns a list of indeces present in both intersection sets which give the index positions
    of SDs that map to an inversion's left breakpoint and have inverted homologous mates that map to the
    inversion's right breakpoint

    :param SD1bed: BedTool of SDs for set 1 made with `_SDdf_to_BedTools`
    :param SD2bed: BedTool of SDs for set 2 made with `_SDdf_to_BedTools`
    :param brk1bed: BedTool object of breakpoint 1 made with `_brkpnt_bed_spans_from_inv`
    :param brk2bed: BedTool object of breakpoint 2 made with `_brkpnt_bed_spans_from_inv`

    :type SD1bed: pybedtools.bedtool.BedTool
    :type SD2bed: pybedtools.bedtool.BedTool
    :type brk1bed: pybedtools.bedtool.BedTool
    :type brk2bed: pybedtools.bedtool.BedTool

    :returns: list of index values of rows in the SD dataframe made with `_filter_SD_bed` in which an SD in its inverted homolog mate map to an inversion's breakpoints
    :rtype: list
    """

    #finding SD intersections with breakpoint 1
    SD_to_brk1 = SD1bed.intersect(brk1bed, loj=True)

    #finding SD intersections with breakpoint 2
    SD_to_brk2 = SD2bed.intersect(brk2bed, loj=True)

    #find idxs for both sets
    sd1_idxs = _find_bedintersect_idxs(SD_to_brk1)
    sd2_idxs = _find_bedintersect_idxs(SD_to_brk2)

    #find idx values present in both sets (intersection)
    sd_intersects_idxs = _find_idxs_in_both_sets(sd1_idxs, sd2_idxs)

    #delete temporary bedfiles
    _delete_BedTool(SD_to_brk1)
    _delete_BedTool(SD_to_brk2)

    return sd_intersects_idxs



def _count_span_results(sddf, intersectidxs):
    """
    Helper function that takes in the filtered SD dataframe and the list of indeces for
    rows of the dataframe where the inversion breakpoints map to inverted SD mates and 
    returns 1) the number of inverted SD pairs and 2) the total number of inverted homologous
    base pairs

    :param sddf: dataframe of filtered bedfile from `_filter_SD_bed`
    :param intersectidxs: list of index values of rows in the SD dataframe from `_find_brkpnt_intersect_indeces`

    :type sddf: pandas.DataFrame
    :type intersectidxs: list

    :returns: number of inverted SD pairs, number of inverted homologous bps
    :rtype: int, int
    """

    #counting number of inverted SD pairs
    num_inv_SDs = len(intersectidxs)

    #counting number of inverted homologus bps
    if num_inv_SDs == 0:
        num_inv_bps = 0
    else:
        filtered_df = sddf.iloc[intersectidxs]
        num_inv_bps = int(filtered_df["matchB"].sum())

    return num_inv_SDs, num_inv_bps



def _randomize_inversion_spans(invbed, chrom_sizes_file, seed_value):
    """
    Helper function take in a BedTool object of inversion spans, randomly shuffles the
    coordinates within the confines of hg38 genome and preserving chromosome and returns 
    each inversion in a list of tuples where each tuple contains the chrom, start, and end of each span

    :param invbed: BedTool object of inversion spans
    :param chrom_sizes_file: file name for tab separated table of chromosome sizes
    :param seed_value: seed value to be used by pybedtools.BedTool.shuffle

    :type invbed: pybedtools.bedtool.BedTool
    :type chrom_sizes_file: str
    :type seed_value: int

    :returns: list of tuples where each tuple is an inversion and contains (chrom, start, end)
    :rtype: [tuple]
    """

    #randomize spans
    randomized_bed = invbed.shuffle(g=chrom_sizes_file, chrom=True, seed=seed_value)

    #convert BedTool to list of tuples
    randomized_list = []
    for interval in randomized_bed:
        random_tuple = (interval[0], int(interval[1]), int(interval[2]))
        randomized_list.append(random_tuple)

    #delete temporary bedfile
    _delete_BedTool(randomized_bed)

    return randomized_list



def _compute_pvals(test_value, null_set):
    """
    Computes a p-value for the test value from a poisson distribution with a lamda = mean(null set),
    and computes an empirical p-value for the test value but finding what proportion of the null set
    is at least as extreme as the test value (one-tailed)

    :param top_value: integer or float for the average value from the set of interest
    :param null_set: list of integers or floats that of average values in the null set from the permutation

    :type top_value: float
    :type null_set: list of floats
    
    :returns: poisson_p, gamma_p, and empirical_p
    :rtype: float, float, float
    """

    #if data is >=1 use poisson, if not use gamma
    if test_value >= 1:
        #estimate lambda
        l = mean(null_set)
        #estimate poisson p-value
        pois_p = scipy.stats.poisson.sf(test_value-1, l)
        #no output for gamma
        gamma_p = np.nan 
    else:
        #fit gamma distribution
        a, loc, scale = scipy.stats.gamma.fit(null_set)
        #estimate gamma p-value
        gamma_p = scipy.stats.gamma.sf(test_value, a, loc, scale)
        #no output for poisson
        pois_p = np.nan

    #compute empirical p-value
    total = len(null_set)
    counter = 0
    for val in null_set:
        if val >= test_value:
            counter += 1
    emp_p = (counter+1)/(total+1)

    return pois_p, gamma_p, emp_p




def _count_results_for_set(sddf, SD1bed, SD2bed, span_set, either_side_bp):
    """
    Function that takes in the filtered SD dataframe from `_filter_SD_bed` and a list of spans
    (either for inversions or random regions) and iterates through them. For each span the function
    finds SD intersects with the breakponts of each span and counts how many inverted SDs intersect
    the breakpoints and how many inverted homologous bps this represents. Then the function finds the 
    average number of intersecting inverted SDs and the average number of inverted homologous bps for
    the set of spans.

    :param sddf: dataframe of filtered bedfile from `_filter_SD_bed`
    :param SD1bed: BedTool of SDs for set 1 made with `_SDdf_to_BedTools`
    :param SD2bed: BedTool of SDs for set 2 made with `_SDdf_to_BedTools`
    :param span_set: list of tuples where each tuple is an inversion and contains (chrom, start, end) from `_load_inversion_bed`
    :param either_side_bp: number of base pairs on either side of each span edge

    :type sddf: pandas.DataFrame
    :type SD1bed: pybedtools.bedtool.BedTool
    :type SD2bed: pybedtools.bedtool.BedTool
    :type span_set: list
    :type either_side_bp: int

    :returns: average number of intersecting inverted SDs in the set, average number of inverted homologous bps in the set
    :rtype: float, float
    """

    #empty lists to fill with span results
    invert_SDs_list = []
    invert_bps_list = []

    #iterate through spans (inv or random)
    for span in span_set:

        #find breakpoint spans
        brk1, brk2 = _brkpnt_bed_spans_from_inv(span, either_side_bp)

        #find row indeces of intersects across span breakpoints
        intersect_row_idxs = _find_brkpnt_intersect_indeces(SD1bed, SD2bed, brk1, brk2)

        #cound results for span
        num_invert_SDs, num_invert_bps = _count_span_results(sddf, intersect_row_idxs)

        #append to lists
        invert_SDs_list.append(num_invert_SDs)
        invert_bps_list.append(num_invert_bps)

        #deleting temporary bedfiles (for each iteration)
        _delete_BedTool(brk1)
        _delete_BedTool(brk2)

    #averaging
    avg_invert_SDs = mean(invert_SDs_list)
    avg_invert_bps = mean(invert_bps_list)


    return avg_invert_SDs, avg_invert_bps



def _delete_BedTool(BedTool_obj):
    """
    Helper function for deleting the temperorary file for a BedTool object

    :param BedTool_obj: BedTool object instance
    :type BedTool_obj: pybedtools.bedtool.BedTool
    """

    fname = BedTool_obj.fn

    if os.path.exists(fname):
        os.remove(fname)





### CLASSES ###

class ArgLoader():
    def __init__(self, SD_bedfile, inversion_bedfile, chrom_sizes_file, num_bp_either_side, seed_val):
        """
        Class used to store values to be passed to the `run_permutation` function as a single
        argument so that this function can be run in parallel with `Pool.map`

        :param SD_bedfile: hg38.chr_only.SDs.bed file path
        :param inversion_bedfile: bedfile name of inversion spans
        :param chrom_sizes_file: file name for tab separated table of chromosome sizes
        :param num_bp_either_side: number of base pairs on either side of each span edge
        :param seed_val: seed value used by `_randomize_inversion_spans`

        :type SD_bedfile: str
        :type inversion_bedfile: str
        :type chrom_sizes_file: str
        :type num_bp_either_side: int
        :type seed_val: int

        Attributes:
            SDbed (str): hg38.chr_only.SDs.bed file path
            Invbed (str): bedfile name of inversion spans
            chrom_lens (str): file name for tab separated table of chromosome sizes
            bp (int): number of base pairs on either side of each span edge
            seed (int): seed value used to randomize inversion spans
        """

        self.SDbed = SD_bedfile
        self.Invbed = inversion_bedfile
        self.chrom_lens = chrom_sizes_file
        self.bp = num_bp_either_side
        self.seed = seed_val



class SetResult():
    def __init__(self, SD_avg, bp_avg):
        """
        Class to store the result of a permutation run or test set run.

        :param SD_avg: average number of inverted SDs counted for a run
        :param bp_avg: average number of inverted homologous bps counted for a run

        :type float:
        :type float:

        Attributes:
            avg_SDs (float): average number of inverted SDs
            avg_bps (float): average number of inverted homologous bps
        """

        self.avg_SDs = SD_avg
        self.avg_bps = bp_avg



class Permutations():
    def __init__(self, list_of_SetResult_objects):
        """
        Class to store the results of each permutation used in `run_the_program`

        :param list_of_SetResult_objects: list of results from each permutation, see `SetResult`

        :type list_of_SetResult_objects: list

        Attributes:
            SD_null_set (list): list of permutation results for # of inverted SDs
            bp_null_set (list): list of permutation results for # of inverted homologous bps
        """

        self.SD_null_set = []
        self.bp_null_set = []

        for result in list_of_SetResult_objects:
            self.SD_null_set.append(result.avg_SDs)
            self.bp_null_set.append(result.avg_bps)



class PermTestResults():
    def __init__(self, test_val, null_set_vals):
        """
        Class to store the results of a permutation test used in `run_the_program`

        :param test_val: average test set value, e.g. average # of inverted SDs at inversion breakpoints
        :param null_set_vals: list of permutation values, e.g. list of average # of inverted SDs at random breakpoints

        :type test_val: float
        :type null_set_vals: list

        Attributes:
            val (float): test value
            null_set (list): list of permutation results
            p_poisson (float): p-value estimated from a poisson distribution
            p_gamma (float): p-value estimated from a gamma distribution
            p_empirical (float): p_value computed as proportion of null values at least as extreme as test value
        """

        #setting val and null_set attributes
        self.val = test_val
        self.null_set = null_set_vals

        #computing Z-score and p-value
        self.p_poisson, self.p_gamma, self.p_empirical = _compute_pvals(self.val, self.null_set)





### FUNCTIONS THAT RUN COMPUTATIONS ON INVERSION SET AND ON PERMUTATONS ###

def inversion_results(SD_bedfile, inversion_bedfile, num_bp_either_side):
    """
    Function computes the average number of inverted SDs and average number of inverted homologous 
    base pairs that map to inversion breakpoints.

    :param SD_bedfile: hg38.chr_only.SDs.bed file path
    :param inversion_bedfile: bedfile name of inversion spans
    :param num_bp_either_side: number of base pairs on either side of each span edge

    :type SD_bedfile: str
    :type inversion_bedfile: str
    :type num_bp_either_side: int

    :returns: Avg # of inverted SDs at inversion breakpoints and Avg # of inverted homologous bps at inversion breakpoints stored as a SetResult object
    :rtype: SetResult
    """

    #loading bedfile of SDs into a filtered dataframe, see `_filter_SD_bed`
    beddf = _filter_SD_bed(SD_bedfile)
    #convert SD1 set and SD2 set to BedTool objects
    sd1, sd2 = _SDdf_to_BedTools(beddf)
    #loading bedfile of inversion spans into a list of tuples, see `_load_inversion_bed`
    invs_list = _load_inversion_bed(inversion_bedfile)
    #compute results for inversion set
    inv_avg_SDs, inv_avg_bps = _count_results_for_set(beddf, sd1, sd2, invs_list, num_bp_either_side)

    inv_result = SetResult(inv_avg_SDs, inv_avg_bps)

    #delete temporary bedfiles
    _delete_BedTool(sd1)
    _delete_BedTool(sd2)

    return inv_result




def run_permutation(argloader_object):
    """
    Function computes the average number of inverted SDs and average number of inverted homologous 
    base pairs that map to inversion breakpoints for a permutation set generated by randomly shuffling
    the BedTool of inversion spans.

    :param argloader_object: object that stores the hg38.chr_only.SDs.bed filename, inversion spans bedfile name, chrom lengths filename, number of base pairs on either side of each span edge, and a seed value

    :type argloader_object: ArgLoader

    :returns: Avg # of inverted SDs and Avg # of inverted homologous bps stored as a SetResult object
    :rtype: SetResult
    """

    SD_bedfile = argloader_object.SDbed
    inversion_bedfile = argloader_object.Invbed
    chromlens = argloader_object.chrom_lens
    num_bp_either_side = argloader_object.bp
    seed = argloader_object.seed

    #loading bedfile of SDs into a filtered dataframe, see `_filter_SD_bed`
    beddf = _filter_SD_bed(SD_bedfile)
    #convert SD1 set and SD2 set to BedTool objects
    sd1, sd2 = _SDdf_to_BedTools(beddf)

    #loading inversions as a BedTool object
    invs_bed = pybedtools.BedTool(inversion_bedfile)

    #randomize set of inversions
    randomized_spans = _randomize_inversion_spans(invs_bed, chromlens, seed)
    #use beddf, sd1, sd2, and the ranomdized set and compute results
    randomized_avg_SDs, randomized_avg_bps = _count_results_for_set(beddf, sd1, sd2, randomized_spans, num_bp_either_side)

    randomized_result = SetResult(randomized_avg_SDs, randomized_avg_bps)

    #delete temporary bedfiles
    _delete_BedTool(sd1)
    _delete_BedTool(sd2)

    return randomized_result





### PLOTTING ###

def plot_permutation_test(perm_result_obj, param_name, plot_title):
    """
    Plots the null distribution in a histogram and the value of interest as a vertical line. Saves a png of plot
    
    :param perm_result_obj: permutation test results, instance of PermTestResults
    :param param_name: name of parameter being plotted, used for x-axis label
    :param plot_title: title of plot, also used as the filename for output png

    :type perm_result_obj: PermTestResults
    :type param_name: str
    :type plot_title: str
    """

    #number of bins using Freedman–Diaconis rule
    q25, q75 = np.percentile(perm_result_obj.null_set, [25, 75])
    bin_width = 2 * (q75 - q25) * len(perm_result_obj.null_set) ** (-1/3)
    if bin_width == 0:
        bins = 10
    else:
        try:
            bins = round((max(perm_result_obj.null_set) - min(perm_result_obj.null_set)) / bin_width)
        except OverflowError:
            bins = 10

    #plotting
    plt.figure(figsize=(8,6))
    plt.hist(perm_result_obj.null_set, bins=bins) #density=True
    plt.ylabel('Density')
    plt.xlabel(param_name)
    plt.title(plot_title)
    plt.axvline(x=perm_result_obj.val, color='r', lw=3)
    plt.text(0.5, 0.8,
             f"Poisson p-value = {round(perm_result_obj.p_poisson, 5)}\nGamma p-value = {round(perm_result_obj.p_gamma, 5)}\nEmpirical p-value = {round(perm_result_obj.p_empirical, 5)}", 
             horizontalalignment='center', verticalalignment='bottom', bbox=dict(facecolor='yellow', alpha=0.25),
             transform=plt.gca().transAxes)
    plt.savefig(f"{plot_title}.png", bbox_inches='tight')
    plt.clf()





### RUNNING THE PROGRAM ###

def run_the_program(SD_bed, inv_bed, chrom_sizes, num_bp, num_perms, num_procs):
    """
    Main function that runs a permutation test between inversions and random regions
    counting the number of inverted SDs and inverted homologous bps between the two
    sets and plots the results.

    :param SD_bed: hg38.chr_only.SDs.bed file path
    :param inv_bed: bedfile name of inversion spans
    :param chrom_sizes: filename of tab separated Chroms and Lengths e.g. chr1    248956422
    :param num_bp: number of base pairs on either side of each span edge
    :param num_perms: number of random permutations of random spans to run
    :param num_procs: number of processes to run permutations in parallel

    :type SD_bed: str
    :type inv_bed: str
    :type chrom_sizes: str
    :type num_bp: int
    :type num_perms: int
    :type num_procs: int
    """

    #Results for inversion set
    inversion_result = inversion_results(SD_bed, inv_bed, num_bp)

    #Generate Arguments for parallel permutation calls
    argloader_list = []
    for i in range(num_perms):
        random_seed = _generate_seed_value()
        argloader_obj = ArgLoader(SD_bed, inv_bed, chrom_sizes, num_bp, random_seed)
        argloader_list.append(argloader_obj)

    #run permutations in parallel
    with Pool(processes=num_procs) as pool:
        null_results = pool.map(run_permutation, argloader_list)


    #collecting results into a Permutations object
    permutation_results = Permutations(null_results)
    
    #Permutation Test Results
    sdpermresults = PermTestResults(inversion_result.avg_SDs, permutation_results.SD_null_set)
    bppermresults = PermTestResults(inversion_result.avg_bps, permutation_results.bp_null_set)

    #ploting results
    plot_permutation_test(sdpermresults, "Avg. Number of Inverted SDs", "Avg. Number of Inverted SD Homologs\nat Inversion Breakpoints vs. Random Regions")
    plot_permutation_test(bppermresults, "Avg. Number of Inverted Homologous BPs", "Avg. Number of Inverted Homologous BPs\nat Inversion Breakpoints vs. Random Regions")

    #cleanup any temp files leftover
    pybedtools.cleanup(remove_all=True)





def main():
    run_the_program("hg38.chr_only.SDs.bed", "porubsky_inversions.bed", "hg38_chrom_sizes.txt", 10000, 10000, 10)

if __name__ == "__main__":
    main()

