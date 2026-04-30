"""
bootstrap_LD v1.1

Dependencies:
- PLINK v1.9
- bcftools 1.23
- grep

Estimates linkage disequilibrium (|D'|) between sites in one vcf.gz file compared to sites in another vcf.gz file
This version does not actually run sampling bootstraps but for now uses all samples in the vcf
The two vcfs must include the same samples.

This new version now reports LD as 1) mean of allm site comparisons, 2) mean of all non-zero site comparisions and 3)geometric mean of all non-zero site comparisions
"""

import subprocess
import tempfile
import gzip
import os
import pandas as pd
import random
from statistics import geometric_mean, mean
random.seed(754)

def _pos2idx(listuniqpos):
    """
    Helper function takes in a list of unique position values and returns a dictionary of {pos:idx}

    :param listuniqpos: List of unique position values from a plink LD output table
    :type listuniqpos: list

    :returns: dictionary of position values mapped to integer (index) values
    :rtype: dict
    """

    #creating dictionary
    pos2idx = {}
    for idx, pos in enumerate(listuniqpos):
        pos2idx[pos] = idx

    return pos2idx


def _open_LDfile_as_df(LDfile):
    """
    Helper function loads in plink LD output and returns the table as a dataframe.

    :param LDfile: plink output LD table file name
    :type LDfile: str

    :returns: dataframe
    :rtype: pandas.DataFrame
    """

    #loading file
    df = pd.read_csv(LDfile, sep=r'\s+')
    df.rename(columns={"DP": "Dprime"}, inplace=True)

    #selecting columns
    df = df[["BP_A", "BP_B", "Dprime"]]

    return df


def _convert_pos_to_idx(LDfile):
    """
    Helper function loads in plink LD output and prepares for plotting by converting position values
    to integer (index) values. Note that all SNPs compared need to be from the same chromosome.

    :param LDfile: plink output LD table file name
    :type LDfile: str

    :returns: dataframe with the converted position values
    :rtype: pandas.DataFrame
    """
    
    #loading file
    df = _open_LDfile_as_df(LDfile)

    #creating position to idx map
    unique_positions_A = list(set(df["BP_A"].to_list()))
    unique_positions_B = list(set(df["BP_B"].to_list()))
    unique_positions_A.sort()
    unique_positions_B.sort()
    pos2idx_dict_A = _pos2idx(unique_positions_A)
    pos2idx_dict_B = _pos2idx(unique_positions_B)
    
    #adding idx_columns
    df["idx_A"] = df["BP_A"].map(pos2idx_dict_A)
    df["idx_B"] = df["BP_B"].map(pos2idx_dict_B)

    #grabbing necessary columns only
    df = df[["idx_A", "idx_B", "Dprime"]]
    
    return df



def _geometric_mean_of_Dprime(dprime_col):
    """
    Helper function takes in a dataframe column of Dprime values and computes the
    geometric mean after clipping out an 0 values.

    :param dprime_col: Dataframe column of Dprime values
    :type dprime_col: pandas.Series

    :returns: Geometric Mean of Dprime values
    :rtype: float
    """

    #clipping out zeros
    dprimes = [x for x in dprime_col if x > 0]

    #computign geometric mean
    gmean_dprime = geometric_mean(dprimes)

    return gmean_dprime



def _means_of_Dprime(dprime_col):
    """
    Helper function takes in a dataframe column of Dprime values and computes the
    mean  and mean after clipping out an 0 values.

    :param dprime_col: Dataframe column of Dprime values
    :type dprime_col: pandas.Series

    :returns: Mean of Dprime values,  Mean of non-zero Dprime values, number of site comparisions, number of non-zero site comparisions 
    :rtype: float, float, int, int
    """


    #clipping out zeros
    dprimes = dprime_col.to_list()
    dprimes_nonzero = [x for x in dprime_col if x > 0]

    #counting site comparisons
    sitecomps = len(dprimes)
    sitecomps_nonzero = len(dprimes_nonzero)

    #computign mean
    mean_dprime = mean(dprimes)
    mean_nonzero_dprime = mean(dprimes_nonzero)

    return mean_dprime, mean_nonzero_dprime, sitecomps, sitecomps_nonzero




def run_LD(vcfgz1, vcfgz2, save_for_plotting=False):
    """
    Function computes LD (|D'|) between vcf1 and vcf2, computes the geometric mean
    of Dprime, and optionally outputs the LD result for plotting (CAUTION THESE FILES CAN BE LARGE)

    :param vcfgz1: filename for first vcf.gz file
    :param vcfgz2: filename for second vcf.gz file
    :param save_for_plotting: boolean value for whether you want to save teh LD table output

    :type vcfgz1: str
    :type vcfgz2: str
    :type save_for_plotting: bool = False

    :returns: geometric mean of Dprime, mean of drpime, mean of non-zero dprime, number of site comparisons, number of non-zero site comparisons, and the LD output name
    :rtype: float, float, float, int, int, str
    """

    print("===================================================================================================================================================")
    print(f"ESTIMATING GEOMETRIC MEAN OF |D'|\nBETWEEN {vcfgz1} AND {vcfgz2}")
    print("===================================================================================================================================================")

    #output file name
    outname = f"{vcfgz1.replace(".vcf.gz", "")}__TO__{vcfgz2.replace(".vcf.gz", ".ld")}"

    #opening files and temp directory
    with tempfile.TemporaryDirectory() as tmpdir:#, open(outname, "w") as outputLD:

        #Temporary file paths
        vcfannA_path = os.path.join(tmpdir, "vcf_ann_A.vcf.gz")
        vcfannB_path = os.path.join(tmpdir, "vcf_ann_B.vcf.gz")
        concatvcf_path = os.path.join(tmpdir, "concat_vcf.vcf.gz")
        idsetA_path = os.path.join(tmpdir, "id_set_A.txt")
        # pruned_out_base = os.path.join(tmpdir, "pruned")
        # pruned_out_path = pruned_out_base + ".prune.out"
        variant_ids_path = os.path.join(tmpdir, "variants.txt")
        trimmed_out_variants_path = os.path.join(tmpdir, "excluded_variants.txt")
        unfilteredLD_path_base = os.path.join(tmpdir, "unfilted_LD")
        unfilteredLD_path = unfilteredLD_path_base + ".ld"
        temp_output = os.path.join(tmpdir, f"TEMP_{outname}")


        #1) Set site IDs
        print("\nAdding SNP IDs to variants in input vcfs...\n")
        subprocess.Popen(["bcftools", "annotate", "--set-id", "SetA:%CHROM:%POS", vcfgz1, "-Oz", "-o", vcfannA_path]).wait()
        subprocess.Popen(["bcftools", "annotate", "--set-id", "SetB:%CHROM:%POS", vcfgz2, "-Oz", "-o", vcfannB_path]).wait()
        subprocess.Popen(["tabix", "-p", "vcf", vcfannA_path]).wait()
        subprocess.Popen(["tabix", "-p", "vcf", vcfannB_path]).wait()

        #2) Query "SetA:%CHROM%POS" > txt
        print("\nGenerating list of SNP IDs in SetA...\n")
        with open(idsetA_path, "w") as setAfile:
            subprocess.Popen(["bcftools", "query", "-f", "SetA:%CHROM:%POS", vcfgz1], stdout=setAfile).wait()

        #3) Concatenate vcfs
        print("\nConcatenating set vcfs for plink input...\n")
        subprocess.Popen(["bcftools", "concat", vcfannA_path, vcfannB_path, "-Oz", "-o", concatvcf_path]).wait()
        subprocess.Popen(["tabix", "-p", "vcf", concatvcf_path]).wait()
        with open(variant_ids_path, "w") as variant_ids:
            subprocess.Popen(["bcftools", "query", "-f", "%ID", concatvcf_path], stdout=variant_ids).wait()

        # #4) Prune Variants
        # print("\nPruning variants with r^2 > 0.5 in 10Kb windows...\n")
        # subprocess.Popen(["plink", "--vcf", concatvcf_path, "--indep-pairwise", "10kb", "1", "0.5", "--out", pruned_out_base]).wait()

        #4) Trimming Variants
        downsample_size = 1000
        print(f"\nDownsampling to {downsample_size} SNPs...\n")
        with open(variant_ids_path, "r") as variants, open(trimmed_out_variants_path, "w") as excluded_variants:
            variants_list = variants.readlines()
            Avariants = [x for x in variants_list if "SetA" in x]
            Bvariants = [x for x in variants_list if "SetB" in x]
            excluded_Avariants_list = random.sample(Avariants, max((len(Avariants)-downsample_size), 0))
            excluded_Bvariants_list = random.sample(Bvariants, max((len(Bvariants)-downsample_size), 0))
            excluded_variants_list = excluded_Avariants_list + excluded_Bvariants_list
            excluded_variants.write("".join(excluded_variants_list))

        #5) Run plink LD with SetA site list
        print("\nEstimating LD between SetA SNPs and all SNPs in concatenated vcf...\n")
        # subprocess.Popen([
        #     "plink", "--vcf", concatvcf_path, 
        #     "--r2", "dprime", "--exclude", pruned_out_path, "--ld-snp-list", idsetA_path,
        #     "--ld-window-kb", "1000000", "--ld-window", "1000000", "--ld-window-r2", "0", 
        #     "--out", unfilteredLD_path_base
        #     ]).wait()
        subprocess.Popen([
            "plink", "--vcf", concatvcf_path, 
            "--r2", "dprime", "--exclude", trimmed_out_variants_path, "--ld-snp-list", idsetA_path,
            "--ld-window-kb", "1000000", "--ld-window", "1000000", "--ld-window-r2", "0", 
            "--out", unfilteredLD_path_base
            ]).wait()

        #6) Grep for SetB:
        print("\nFiltering for SetA to SetB comparisons...\n")
        with open(temp_output, "w") as tempout:
            subprocess.Popen(["grep", "-E", "CHR_A|SetB:", unfilteredLD_path], stdout=tempout).wait()

        #7) Compute Geometric Mean of D'
        if save_for_plotting is True:
            #7a) Convert position coordinates to indeces and save for plotting
            print("\nConverting position coordinates to index values...\n")
            converted_LD_df = _convert_pos_to_idx(temp_output)
            converted_LD_df.to_csv(outname, index=False)
            print("\nComputing Geometric Mean of Dprime...\n")
            gm = _geometric_mean_of_Dprime(converted_LD_df["Dprime"])
            m, m_nz, num_site_comps, num_site_comps_nz = _means_of_Dprime(LD_df["Dprime"])
        else:
            LD_df = _open_LDfile_as_df(temp_output)
            print("\nComputing Geometric Mean of Dprime...\n")
            gm = _geometric_mean_of_Dprime(LD_df["Dprime"])
            m, m_nz, num_site_comps, num_site_comps_nz = _means_of_Dprime(LD_df["Dprime"])

    return gm, m, m_nz, num_site_comps, num_site_comps_nz, outname





def _parse_comparison_file(comparison_file):
    """
    Helper function that parses the comparison file and returns a list of comparisons.
    The list will be a list of lists with each sublist being a comparison to input into
    `run_LD`. Comparisons should look like ["file1.vcf.gz", "file2.vcf.gz"]
    Comparisons that are to be plotted should look like ["file1.vcf.gz", "file2.vcf.gz", "PLOT"]

    :param comparisons_file: name of file of list of comparisons to run
    :type comparisons_file: str

    :returns: list of comparisons
    :rtype: [list]
    """

    #comparisons
    comparison_list = []

    #opening file
    with open(comparison_file, "r") as cfile:
        for line in cfile:
            line_list = line.strip().split(",")
            comparison_list.append(line_list)

    return comparison_list


def _parse_split_name(name):
    """
    Helper function used by `_parse_outname` that takes a peak name i.e. "peak30.chr15_29990002_30910000.invbrk6.chr15_30618104_NA"
    and parses it to find the peak name, coordinate, and feature type

    :param name: name from LD output filename before or after the "__TO__"
    :type name: str

    :returns: peak name, chr, start, end, and feature type
    :rtype: str, str, int, int, str
    """

    #split name
    split_peak_name = name.split(".")

    #grab data
    peakname = split_peak_name[0]
    coord_str = split_peak_name[1]
    featurename = split_peak_name[2]

    #clean data
    split_coord = coord_str.split("_")
    chrom = split_coord[0]
    start = int(split_coord[1])
    end = int(split_coord[2])

    if featurename == "invbrk0":
        featurename = "None"
    else:
        featurename = featurename.replace("brk", "")

    return peakname, chrom, start, end, featurename




def _parse_outname(outname):
    """
    Helper function used by `per_comparison_dprime_gmeans` that parses the LD output file name (whether or not this is actually being saved)
    to find the name, coordinate, and feature type of peakA and peakB to be put into the output table

    :param outname: LD output name
    :type outname: str

    :returns: dictionary of name, coordinates, and feature type of peakA and peakB, comparison name and comparison_type
    :rtype: {"peakA":str, "chrA":str, "startA":int, "endA":int, "featureA":str, "peakB":str, "chrB":str, "startB":int, "endB":int, "featureB":str, "comparison_name":str, "comparison_type":str}
    """

    #initialize empty dictionary
    parsed_dict = {}

    #split name into A and B
    split_name = outname.split("__TO__")
    nameA = split_name[0]
    nameB = split_name[1]

    #parse A and B
    peak_A, chr_A, start_A, end_A, feature_A = _parse_split_name(nameA)
    peak_B, chr_B, start_B, end_B, feature_B = _parse_split_name(nameB)

    #comparison name
    comp_name = f"{peak_A} to {peak_B}"

    #finding comparison type
    if feature_A != "None" and feature_A == feature_B:
        comp_type = "INV to INV"
    elif "inv" in feature_A and "inv" in feature_B:
        comp_type = "Adjacent INVs"
    elif feature_A == "None" and feature_B == "None":
        comp_type = "None to None"
    else:
        comp_type = "INV to None"

    #add to dictionary
    parsed_dict["peakA"] = peak_A
    parsed_dict["chrA"] = chr_A
    parsed_dict["startA"] = start_A
    parsed_dict["endA"] = end_A
    parsed_dict["featureA"] = feature_A
    parsed_dict["peakB"] = peak_B
    parsed_dict["chrB"] = chr_B
    parsed_dict["startB"] = start_B
    parsed_dict["endB"] = end_B
    parsed_dict["featureB"] = feature_B
    parsed_dict["comparison_name"] = comp_name
    parsed_dict["comparison_type"] = comp_type

    return parsed_dict



def per_comparison_dprime_gmeans(comparisons_file, outcsv):
    """
    Function to run `run_LD` on a list of comparisons and generate a table of geometric mean
    for each comparision. The function will also save the LD output for specified comparisons 
    you want to plot.

    Note on comparisons_file:
        - each line is a comparison
        - items in each line must be comma-separated i.e. "file1.vcf.gz,file2.vcf.gz"
        - if you would like a comparison to be plotted add comma + "PLOT" as the end of the line i.e. "file1.vcf.gz,file2.vcf.gz,PLOT"
        - This function will check if a line contains "PLOT" and 3 elements to decide if output LD should be saved

    :param comparisons_file: name of file of list of comparisons to run
    :param outcsv: base name of output file

    :type comparisons_file: str
    :type outcsv: str
    """

    #loading comparisons
    comp_list = _parse_comparison_file(comparisons_file)

    #for outputing table
    label_cols = ["peakA", "chrA", "startA", "endA", "featureA", "peakB", "chrB", "startB", "endB", "featureB", "comparison_name", "comparison_type"]
    outdict = {
        "peakA":[], "chrA":[], "startA":[], "endA":[], "featureA":[], 
        "peakB":[], "chrB":[], "startB":[], "endB":[], "featureB":[], 
        "comparison_name":[], "comparison_type":[], "Num_Site_Comparisons":[], "Num_Site_Comparisons_nonzero":[],
        "dprime_gmean":[], "dprime_mean":[], "dprime_mean_nonzero":[]
        }

    #deciding what to do with each comparison
    for comp in comp_list:

        #no plotting
        if len(comp) == 2:
            geometricmean, avg, nonzeroavg, num_site_comparisons, num_site_comparisons_nonzero, comp_label = run_LD(comp[0], comp[1])
            parsed_outname = _parse_outname(comp_label)
            for col in label_cols:
                outdict[col].append(parsed_outname[col])
            outdict["dprime_gmean"].append(geometricmean)
            outdict["dprime_mean"].append(avg)
            outdict["dprime_mean_nonzero"].append(nonzeroavg)
            outdict["Num_Site_Comparisons"].append(num_site_comparisons)
            outdict["Num_Site_Comparisons_nonzero"].append(num_site_comparisons_nonzero)

        #yes plotting
        elif len(comp) == 3 and comp[2] == "PLOT":
            geometricmean, avg, nonzeroavg, num_site_comparisons, num_site_comparisons_nonzero, comp_label = run_LD(comp[0], comp[1], save_for_plotting=True)
            parsed_outname = _parse_outname(comp_label)
            for col in label_cols:
                outdict[col].append(parsed_outname[col])
            outdict["dprime_gmean"].append(geometricmean)
            outdict["dprime_mean"].append(avg)
            outdict["dprime_mean_nonzero"].append(nonzeroavg)
            outdict["Num_Site_Comparisons"].append(num_site_comparisons)
            outdict["Num_Site_Comparisons_nonzero"].append(num_site_comparisons_nonzero)


        else:
            continue

    #contructing table
    df = pd.DataFrame(outdict)
    df.to_csv(f"{outcsv}.csv", index=False)






per_comparison_dprime_gmeans("comparisons.txt", "means_dprime_comparisons_with_all_samples")

