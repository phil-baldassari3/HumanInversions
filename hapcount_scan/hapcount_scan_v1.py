"""
This script takes in a phased and imputed vcf of biallelic SNPs for a single chromosome
and counts how many unique haplotype are contained within each window of defined number
of SNPs that slides for a defined number of SNPs. The program outputs a bedfile with
the windowed results. This can be run on multiple vcfs at once using mutliprocessing.
"""

import gzip
import pandas as pd
import numpy as np
from collections import deque
from multiprocessing import Pool


### HELPER FUNCTIONS ###

def _joinany(sep, mylist):
    """
    Helper function to join a list by a sep regardless of element types
    """

    strs = [str(x) for x in mylist]

    return sep.join(strs)



def _pop_step(d, step):
    """
    Helper function used to pop a number of elements from the left
    of a deque equal to the step size in. This is done in place.
    """

    for i in range(step):
        d.popleft()


def _construct_matrix(d):
    """
    Helper function takes in the window deque and constructs the haplo matrix
    
    :param d: window deque
    :type d: deque of RecordLoader objects

    :returns: each sublist is a site and each sublist element is the haplotype at that site
    :rtype: list of lists
    """
 

    matrix = []
    for record in d:
        hap_ls = []
        for gt in record.GTs:
            haps = gt.split("|")
            for h in haps:
                hap_ls.append(int(h))
        matrix.append(hap_ls)


    return matrix



def _count_unique_haps(haplo_matrix):
    """
    Helper counts unique haplotypes and sample size of a given window.
    
    :param haplo_matrix: each sublist is a site and each sublist element is the haplotype at that site
    :type haplo_matrix: list of lists

    :returns: number of unique haplotypes and total number of haplotypes (sample size) in the haplo_matrix
    :rtype: int, int
    """


    #converting matrix to array and transposing
    matrix = np.array(haplo_matrix)
    matrix = np.transpose(matrix)

    #counting sample size
    num_samples = len(matrix)

    #uniquify
    uniq_matrix = np.unique(matrix, axis=0)

    #counting unique haplotypes
    num_uniq_haps = len(uniq_matrix)


    return num_uniq_haps, num_samples



### CLASSES ###

class RecordLoader():
    def __init__(self, vcfline):
        """
        This class is for parsing a vcf record for it to be added to a deque.
        
        :param vcfline: record line from vcf to be parsed
        :type vcfline: str

        Attributes:
            position (int): record position
            GTs (list of str): list of individual genotypes e.g. ["0|0", "1|0"...]
        """


        self.position = int(vcfline.strip().split("\t")[1])
        self.GTs = vcfline.strip().split("\t")[9:]



class ArgLoader():
    def __init__(self, vcfgz_file, windowsize, windowstep):
        """
        Class used to store arguments for the SNPwindow_hap_counter function which is run using the run_hapcount_scan main function.
        This is done so that the main function takes a single argument so that it can be run in parallel using multiprocessing.Pool
        
        :param vcfgz_file: filename or path to file for the vcf to scan
        :type vcfgz_file: str
        :param windowsize: SNP window size
        :type windowsize: int
        :param windowstep: SNP window step, recommended to be 10% of windowsize
        :type windowstep: int

        Attributes:
            vcf (str): holds the vcf filename for run_hapcount_scan to parse
            winsize (int): holds the window size for run_hapcount_scan to parse
            winstep (int): holds the window step for run_hapcount_scan to parse
        """

        self.vcf = vcfgz_file
        self.winsize = windowsize
        self.winstep = windowstep

    def __str__(self):
        return self.vcf + "\t" + str(self.winsize) + "\t" + str(self.winstep)




### SLIDING WINDOW SCAN FUNCTION ###

def SNPwindow_hap_counter(vcfgz, SNPwindow_size, SNPwindow_step):
    """
    Function counts unique haplotypes in sliding windows on a gizped vcf. See other functions for details.
    
    :param vcfgz: gziped vcf file name or path to file
    :type vcfgz: str
    :param SNPwindow_size: number of SNPs defining the window size
    :type SNPwindow_size: int
    :param SNPwindow_step: number of SNPs to slide window, recommended to be 10% of SNPwindow_size
    :type SNPwindow_step: int
    
    :returns: appends to output bedfile
    :rtype: None
    """

    #setting up deque for sliding window
    window_deque = deque()

    #open files
    chrom = vcfgz.split("_")[-1].replace(".vcf.gz", "")
    with gzip.open(vcfgz, "rt") as vcf, open(f"{chrom}_SNPwindow{SNPwindow_size}_SNPstep{SNPwindow_step}_hap_counts.bed", "a") as outbed:

        #write header line
        outbed.write("\t".join(["#CHROM", "START", "END", "SNP_density", "hapcount", "sample_size"]) + "\n")

        #loop through lines
        for line in vcf:

            #skip header lines
            if line.startswith("#"):
                continue

            #sliding window scan
            else:

                #load record into deque
                window_deque.append(RecordLoader(line))

                #check if window size was reached, if yes, move the window by a step, output results, and continue
                if len(window_deque) >= SNPwindow_size:

                    #grab window position
                    start_pos = window_deque[0].position
                    end_pos = window_deque[-1].position + 1
                    
                    #construct haplo matrix
                    hap_m = _construct_matrix(window_deque)
                    #count haps
                    num_haps, sample_size = _count_unique_haps(hap_m)
                    #compute SNP density
                    snpden = SNPwindow_size / (end_pos - start_pos)
                    #sliding to next step
                    _pop_step(window_deque, SNPwindow_step)

                    #write to file
                    outbed.write(_joinany("\t", [chrom, start_pos, end_pos, snpden, num_haps, sample_size]) + "\n")

                else:
                    continue

        #clean up last window
        if len(window_deque) > 0:
            
            #grab window position
            start_pos = window_deque[0].position
            end_pos = window_deque[-1].position + 1

            #construct haplo matrix
            hap_m = _construct_matrix(window_deque)
            #count haps
            num_haps, sample_size = _count_unique_haps(hap_m)
            #compute SNP density
            snpden = len(window_deque) / (end_pos - start_pos)
            
            #write to file
            outbed.write(_joinany("\t", [chrom, start_pos, end_pos, snpden, num_haps, sample_size]) + "\n")




### MAIN FUNCTION ###
def run_hapcount_scan(argloader_obj):
    """
    Main Function that runs SNPwindow_hap_counter with a single argument of the class ArgLoader
    """

    SNPwindow_hap_counter(argloader_obj.vcf, argloader_obj.winsize, argloader_obj.winstep)



# args = ArgLoader("test1.vcf.gz", 100, 10)
# run_hapcount_scan(args)



### PARALLELIZATION ###

def main():

    vcflist = [
        "beagle_phased_biallelic_SNPs_1000GP30X_chr1.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr2.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr3.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr4.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr5.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr6.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr7.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr10.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr11.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr12.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr13.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr14.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr16.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr18.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr20.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr21.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr22.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_nonPAR_chrX.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_PAR_chrX.vcf.gz",
        "beagle_imputed_males_biallelic_SNPs_1000GP30X_chrY.vcf.gz"
        ]
    
    loaderlist = []

    for vcffile in vcflist:
        loaderlist.append(ArgLoader(vcffile, 1000, 100))
    for vcffile in vcflist:
        loaderlist.append(ArgLoader(vcffile, 10000, 1000))

    # print(len(loaderlist))
    # for l in loaderlist:
    #     print(l)

    pool = Pool(processes=25)
    pool.map(run_hapcount_scan, loaderlist)


if __name__ == '__main__':
    main()
