"""
bootstrap_LD v3

Dependencies:
- PLINK v1.9
- bcftools 1.23
- grep

Estimates linkage disequilibrium (|D'|) between sites in one vcf.gz file compared to sites in another vcf.gz file
for a specified number of boostraps where each boostrap downsamples to 10 individuals.
The two vcfs must include the same samples.

This version now reports 1) mean of all site comparisons, 2) mean of all non-zero site comparisions and 3) geometric mean of all non-zero site comparisions
"""

import subprocess
import tempfile
import gzip
import os
import pandas as pd
import math
import random
from statistics import geometric_mean, mean
from multiprocessing import Pool

random.seed(301)



############################################################################################
##################################### HELPER FUNCTIONS #####################################
############################################################################################

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
    if len(dprimes) == 0:
        gmean_dprime = math.nan
    else:
        gmean_dprime = geometric_mean(dprimes)

    return gmean_dprime


def _means_of_Dprime(dprime_col):
    """
    Helper function takes in a dataframe column of Dprime values and computes the
    mean  and mean after clipping out an 0 values.

    :param dprime_col: Dataframe column of Dprime values
    :type dprime_col: pandas.Series

    :returns: Mean of Dprime values and Mean of non-zero Dprime values
    :rtype: float, float
    """


    #clipping out zeros
    dprimes = dprime_col.to_list()
    dprimes_nonzero = [x for x in dprime_col if x > 0]

    #computign mean
    mean_dprime = mean(dprimes)
    mean_nonzero_dprime = mean(dprimes_nonzero)

    return mean_dprime, mean_nonzero_dprime



def _parse_vcf_name(vcf_filename):
    """
    Helper function that takes the name of an input vcf and parses it to 
    find the peak name, coordinate, and feature type. Note that the vcf
    needs to be named like the following example: "peak3.chr2_29990_30910.invbrk6.chr2_30618104_NA"
    meaning Peak #3, on chr2:29990-30910 that includes inversion #6's left breakpoint at chr2:30618104
    Use "invbrk0" for no inversion. the number after "invbrk" refers to the inversion and not the breakpoint
    The breakpoint is signified by the following coordinate: chr1_1000_NA for left and chr1_NA_2000 for right.

    :param vcf_filename: name of input vcf file
    :type vcf_filename: str

    :returns: peak name, chr, start, end, and feature type
    :rtype: str, str, int, int, str
    """

    #split name
    split_name = vcf_filename.split(".")

    #grab data
    peakname = split_name[0]
    coord_str = split_name[1]
    featurename = split_name[2]

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



def _find_comparison_type(feature1, feature2):
    """
    Helper function that takes the feature names from `_parse_vcf_name` and names the comparison type
    either "INV to INV", "Adjacent INVs", "INV to None" or "None to None"

    :param feature1: feature label
    :param feature2: feature label

    :type feature1: str
    :type feature2: str

    :returns: comparison type
    :rtype: str
    """


    #finding comparison type
    if feature1 != "None" and feature1 == feature2:
        comp_type = "INV to INV"
    elif "inv" in feature1 and "inv" in feature2:
        comp_type = "Adjacent INVs"
    elif feature1 == "None" and feature2 == "None":
        comp_type = "None to None"
    else:
        comp_type = "INV to None"

    return comp_type



def _grab_sample_list(vcf_file):
    """
    Helper function used by `get_bootstrap_lists` to get a list of samples contained in a vcf

    :param vcf_file: gziped vcf file path
    :type vcf_file: str

    :returns: list of sample names
    :rtype: list
    """

    with gzip.open(vcf_file, "rt") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                sample_list = line.strip().split("\t")[9:]

    return sample_list


def _parse_comparison_file(comparison_file):
    """
    Helper function that parses the comparison file and returns a list of comparisons.
    The list will be a list of tuples with each tuple containing 2 vcf.gz file names.
    The input file should contain a line for each comparison that looks like: "file1.vcf.gz,file2.vcf.gz"

    :param comparisons_file: name of file of list of comparisons to run
    :type comparisons_file: str

    :returns: list of comparisons
    :rtype: [tuple]
    """

    #comparisons
    comparison_list = []

    #opening file
    with open(comparison_file, "r") as cfile:
        for line in cfile:
            line_list = line.strip().split(",")[:2]
            comparison_list.append(tuple(line_list))

    return comparison_list

############################################################################################





############################################################################################
####################################### CLASSES ############################################
############################################################################################

class TempDataPaths():
    def __init__(self, tmpdir_connection):
        """
        This class takes in the temporary directory connection from `tempfile.TemporaryDirectory` and
        saves temporary file paths used in this program as attributes

        Attributes:
            vcfID_A (str): path to file in temporary directory for vcf of SetA SNPs
            vcfID_B (str): path to file in temporary directory for vcf of SetB SNPs
            concatvcf (str): path to file in temporary directory for concatenated vcf with SNP ID annotations
            ID_setA (str): path to file in temporary directory for list of SetA SNP IDs
            variant_IDs (str): path to file in temporary directory for list of SNP IDs
            trimmed_variants (str): path to file in temporary directory for list of SNPs to exclude
        """

        #paths for vcfs with annntated SNP IDs
        self.vcfID_A = os.path.join(tmpdir_connection, "vcf_ann_A.vcf.gz")
        self.vcfID_B = os.path.join(tmpdir_connection, "vcf_ann_B.vcf.gz")

        #path for concatenated vcf
        self.concatvcf = os.path.join(tmpdir_connection, "concat_vcf.vcf.gz")

        #path for file that lists SNP IDs for SetA
        self.ID_setA = os.path.join(tmpdir_connection, "id_set_A.txt")

        #path for file that list all SNP IDs
        self.variant_IDs = os.path.join(tmpdir_connection, "variants.txt")

        #path for file that lists all SNP IDs excluded in variant thinning
        self.trimmed_variants = os.path.join(tmpdir_connection, "excluded_variants.txt")


class TempBootstrapPaths():
    def __init__(self, tmpdir_connection):
        """
        This class takes in the temporary directory connection from `tempfile.TemporaryDirectory` and
        saves temporary file paths used in this program as attributes

        Attributes:
            bstrap_vcf (str): path to file in temporary directory for subset vcf for bootstrap
            unfilteredLD_base (str): path to file ***basename*** (no extension) in temporary directory for unfiltered LD output
            unfilteredLD (str): path to file in temporary directory for for unfiltered LD output
            filteredLD (str): path to file in temporary directory for filtered (SetA to SetB) LD output
        """

        #temporary bootstrap vcf
        self.bstrap_vcf = os.path.join(tmpdir_connection, "bstrap.vcf.gz")

        #unfiltered plink LD output, path with just base name and path with full file extension
        self.unfilteredLD_base = os.path.join(tmpdir_connection, "unfilted_LD")
        self.unfilteredLD = self.unfilteredLD_base + ".ld"

        #path for file of filtered LD from SetA to SetB
        self.filteredLD = os.path.join(tmpdir_connection, f"TEMP_outputLD")


class ArgLoader():
    def __init__(self, vcf1, vcf2, n_bootstraps, n_thin_SNPs):
        """
        Class used to store arguments for `bootstrapLD_comparison` function which is run using `run_bootstrapLD_comparison` main function.
        This is done so that the main function takes a single argument so that it can be run in parallel using multiprocessing.Pool
        
        :param vcf1: filename or path to file for the vcfA
        :param vcf1: filename or path to file for the vcfB
        :param n_bootstraps: number of bootstraps to run
        :param n_thin_SNPs: number of SNPs to keep in each set

        :type vcf1: str
        :type vcf2: str
        :type n_bootstraps: int
        :type n_thin_SNPs: int

        Attributes:
            vcfA (str): holds the vcf filename for vcf1
            vcfB (str): holds the vcf filename for vcf2
            bstraps (int): holds the number of bootstraps
            nthin (int): holds the number of SNPs to thin to in each set
        """

        self.vcfA = vcf1
        self.vcfB = vcf2
        self.bstraps = n_bootstraps
        self.nthin = n_thin_SNPs

############################################################################################




############################################################################################
################################ Functions to run LD #######################################
############################################################################################

def process_vcfs(vcfgz1, vcfgz2, variant_thin, temppaths_obj):
    """
    Function takes in 2 vcf files to be compared with LD and processes them for input
    into the boostrapping and LD computation functions.

    :param vcfgz1: filename for first vcf.gz file
    :param vcfgz2: filename for second vcf.gz file
    :param variant_thin: variants are thinned out by randomly subsetting so on this number of variants remain in each Set
    :param temppaths_obj: object that contains temporary file paths. See `TempDataPaths`

    :type vcfgz1: str
    :type vcfgz2: str
    :type variant_thin: int
    :type temppaths_obj: TempDataPaths
    """

    #Setting SNP IDs
    ## print("Adding SNP IDs to variants in input vcfs...")
    subprocess.Popen(["bcftools", "annotate", "--set-id", "SetA:%CHROM:%POS", vcfgz1, "-Oz", "-o", temppaths_obj.vcfID_A]).wait()
    subprocess.Popen(["bcftools", "annotate", "--set-id", "SetB:%CHROM:%POS", vcfgz2, "-Oz", "-o", temppaths_obj.vcfID_B]).wait()
    subprocess.Popen(["tabix", "-p", "vcf", temppaths_obj.vcfID_A]).wait()
    subprocess.Popen(["tabix", "-p", "vcf", temppaths_obj.vcfID_B]).wait()

    #writing a list of SetA SNP IDs to file
    ## print("Generating list of SNP IDs in SetA...")
    with open(temppaths_obj.ID_setA, "w") as setAfile:
        subprocess.Popen(["bcftools", "query", "-f", "SetA:%CHROM:%POS", vcfgz1], stdout=setAfile).wait()

    #Concatenating vcfs with annoatated SNP IDs
    ## print("Concatenating set vcfs for plink input...")
    subprocess.Popen(["bcftools", "concat", temppaths_obj.vcfID_A, temppaths_obj.vcfID_B, "-Oz", "-o", temppaths_obj.concatvcf], stderr=subprocess.DEVNULL).wait()
    subprocess.Popen(["tabix", "-p", "vcf", temppaths_obj.concatvcf]).wait()

    #getting list of SNP IDs
    with open(temppaths_obj.variant_IDs, "w") as snp_ids:
        subprocess.Popen(["bcftools", "query", "-f", "%ID", temppaths_obj.concatvcf], stdout=snp_ids).wait()

    #Thinning variants
    ## print(f"Downsampling to {variant_thin} SNPs in each Set...")
    with open(temppaths_obj.variant_IDs, "r") as variants, open(temppaths_obj.trimmed_variants, "w") as excluded_variants:
        variants_list = variants.readlines()
        Avariants = [x for x in variants_list if "SetA" in x]
        Bvariants = [x for x in variants_list if "SetB" in x]
        excluded_Avariants_list = random.sample(Avariants, max((len(Avariants)-variant_thin), 0))
        excluded_Bvariants_list = random.sample(Bvariants, max((len(Bvariants)-variant_thin), 0))
        excluded_variants_list = excluded_Avariants_list + excluded_Bvariants_list
        excluded_variants.write("".join(excluded_variants_list))



def get_bootstrap_lists(num_bootstraps, temppaths_obj):
    """
    Function gets a list of samples form the concatenated vcf whose path is saved in the input
    TempPaths instance and returns a list of bootstrap downsamples to be used by bcftools for
    subsetting

    :param num_bootstraps: number of bootstraps to return
    :param temppaths_obj: object that contains temporary file paths. See `TempDataPaths`

    :type num_bootstraps: int
    :type temppaths_obj: TempDataPaths

    :returns: list of sample lists for each bootstrap. Each sample list is a comma separated str as prescribed by bcftools, e.g. "s1,s2,s3"
    :rtype: list
    """

    #get list of samples
    samples = _grab_sample_list(temppaths_obj.concatvcf)

    #get bootstrap list
    bootstap_list = []
    for i in range(num_bootstraps):
        bootstrap = random.sample(samples, 10)
        bootstap_list.append(",".join(bootstrap))

    return bootstap_list


def LD_on_bootstraps(bootstraps, temppaths_obj):
    """
    Function takes a list of bootstraps, downsamples the concatenated vcf, and runs plink LD between 
    SetA and SetB variants while excluded variants that have been trimmed out

    :param bootstraps: list of bootstraps from `get_bootsrap_lists`
    :param temppaths_obj: object that contains temporary file paths. See `TempDataPaths`

    :type bootstraps: list
    :type temppaths_obj: TempDataPaths

    :returns: boostrap numbers, number of site comparisions, Gmean(|D'|), mean(|D'|), mean(|D'|>0)
    :rtype: list, list, list, list, list
    """

    #setting empty lists for bootstap number and Gmean of D'
    bootstrap_nums = []
    num_site_comparisons = []
    gmean_dprimes = []
    mean_dprimes = []
    mean_nz_dprimes = []

    #iterating through bootstraps
    #total_bootstraps = len(bootstraps)
    for idx, bstrap in enumerate(bootstraps):

        #print statement for tracking progress
        #if idx % 500 == 0:
            ## print(f"Running LD on Bootstrap {idx+1} of {total_bootstraps}...")

        #open another temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            bstrap_paths = TempBootstrapPaths(tmpdir)
            
            #subset bootstrap samples
            subprocess.Popen(["bcftools", "view", "-s", bstrap, temppaths_obj.concatvcf, "-Oz", "-o", bstrap_paths.bstrap_vcf]).wait()

            #run plink LD
            subprocess.Popen([
                "plink", "--silent", "--vcf", bstrap_paths.bstrap_vcf, 
                "--r2", "dprime", "--exclude", temppaths_obj.trimmed_variants, "--ld-snp-list", temppaths_obj.ID_setA,
                "--ld-window-kb", "1000000", "--ld-window", "1000000", "--ld-window-r2", "0", 
                "--out", bstrap_paths.unfilteredLD_base
                ]).wait()
            
            #filter LD
            with open(bstrap_paths.filteredLD, "w") as outputLD:
                subprocess.Popen(["grep", "-E", "CHR_A|SetB:", bstrap_paths.unfilteredLD], stdout=outputLD).wait()

            #compute output stats
            LD_df = _open_LDfile_as_df(bstrap_paths.filteredLD)
            num_site_comps = len(LD_df)
            gm = _geometric_mean_of_Dprime(LD_df["Dprime"])
            m, m_nz = _means_of_Dprime(LD_df["Dprime"])

            #append output stats
            bootstrap_nums.append(idx)
            num_site_comparisons.append(num_site_comps)
            gmean_dprimes.append(gm)
            mean_dprimes.append(m)
            mean_nz_dprimes.append(m_nz)


    return bootstrap_nums, num_site_comparisons, gmean_dprimes, mean_dprimes, mean_nz_dprimes

############################################################################################           
        





############################################################################################
################# Main Function to run bootstrap LD for one comparison #####################
############################################################################################

def bootstrapLD_comparison(vcfA, vcfB, bootstraps, variant_thinning):
    """
    Main Function that runs bootstrap LD between 2 gzipped vcf files. See `_parse_vcf_name` for how the vcf files should be named.
    The function will return a dataframe of Geometric_mean(D') values for each bootstrap along with metadata and  a dataframe of
    the samples included in each bootstrap. A csv with teh results is saved to the working directory.   

    :param vcfA: filename for first vcf.gz file
    :param vcfB: filename for second vcf.gz file
    :param bootstraps: number of bootstraps to run
    :param variant_thinning: number of SNPs to keep in each Set after thinning

    :type vcfA: str
    :type vcfB: str
    :type bootstraps: int
    :type variant_thinning: int

    :returns: Dataframe of Gmean(D')s, Sample lists, and metadata
    :rtype: pandas.DataFrame
    """

    ## print("===================================================================================================================================================")
    ## print(f"ESTIMATING BOOTSTRAPPED GEOMETRIC MEANS OF |D'|\nBETWEEN {vcfA} AND {vcfB}")
    ## print("===================================================================================================================================================")

    #get peak and feature labels for comparision
    peakA, chromA, startA, endA, featureA = _parse_vcf_name(vcfA)
    peakB, chromB, startB, endB, featureB = _parse_vcf_name(vcfB)
    comparison_name = f"{peakA} to {peakB}"
    comparison_type = _find_comparison_type(featureA, featureB)


    #opening temporary directory and defining temporary file paths in the TempDataPaths instance
    with tempfile.TemporaryDirectory() as tmpdir:
        temp_paths = TempDataPaths(tmpdir)

        #process vcfs
        process_vcfs(vcfA, vcfB, variant_thinning, temp_paths)

        #sample for bootstraps
        bootstrap_samples = get_bootstrap_lists(bootstraps, temp_paths)

        #run LD on bootstraps
        bstrap_idxs, site_comp_nums, gm_dprimes, m_dprimes, m_nz_dprimes = LD_on_bootstraps(bootstrap_samples, temp_paths)

        #temporary directory is now closed adn deleted

    #constructing output dfs
    bootstrap_samples_output = [x.replace(",", ":") for x in bootstrap_samples]
    dout = {
        "peakA": peakA, "chrA": chromA, "startA": startA, "endA": endA, "featureA": featureA,
        "peakB": peakB, "chrB": chromB, "startB": startB, "endB": endB, "featureB": featureB,
        "comparison_name": comparison_name, "comparison_type": comparison_type,
        "Num_Site_Comparisons": site_comp_nums, "Bootstrap_ID": bstrap_idxs,
        "dprime_gmean": gm_dprimes, "dprime_mean": m_dprimes, "dprime_mean_nonzero": m_nz_dprimes,
        "Samples": bootstrap_samples_output
    }


    out_df = pd.DataFrame(dout)
    out_df.to_csv(f"{comparison_name.replace(" ", "_")}_bootstrapLD.csv", index=False)

    return out_df

############################################################################################




###### RUNNING THE PROGRAM #####

### MAIN FUNCTIONS ###
def run_bootstrapLD_comparison(argloader_obj):
    """
    Main function runs `bootstrapLD_comparison` with a single argument of the class ArgLoader
    """

    bootstrapLD_comparison(argloader_obj.vcfA, argloader_obj.vcfB, argloader_obj.bstraps, argloader_obj.nthin)



def mainfunc(comp_file, num_bstraps, num_thin_SNPs, n_procs):
    """
    Main function that runs `run_bootstrapLD_comparison` in parallel.

    :param comp_file: filename to file that lists comparisions to be run. See `_parse_comparison_file`
    :param num_bstraps: number of bootstraps to run on each comparision
    :param num_thin_SNPs: number of SNPs to thin each set to
    :param n_procs: number of parallel processes

    :type comp_file: str
    :type num_bstraps: int
    :type num_thin_SNPs: int
    :type n_procs: int
    """
    
    #parse comparison file
    comp_list = _parse_comparison_file(comp_file)

    #load arguments into list of ArgLoader objects
    loader_list = []
    for vcf_file1, vcf_file2 in comp_list:
        loader_list.append(ArgLoader(vcf_file1, vcf_file2, num_bstraps, num_thin_SNPs))

    #run in parallel
    pool = Pool(processes=n_procs)
    pool.map(run_bootstrapLD_comparison, loader_list)


if __name__ == '__main__':
    mainfunc("comparisons.txt", 3000, 1000, 10)
    #mainfunc("comparisons.txt", 10, 100, 10)



### You need to run concat_csvs.py to concatenate the resulting csvs ###