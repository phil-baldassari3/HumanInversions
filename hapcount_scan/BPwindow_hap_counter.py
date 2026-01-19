import gzip
import pandas as pd
import numpy as np
from collections import deque
from multiprocessing import Pool


### HELPER FUNCTIONS ###

def _joinany(sep, mylist):
    """Function to join a list by a sep regardless of element types"""

    strs = [str(x) for x in mylist]

    return sep.join(strs)



def _pop_bp_step(d, step):
    """
    Helper function used to pop a number of elements from the left
    of a deque equal to the step size in. This is done in place.
    """

    #count how many elements to pop
    count_ele_to_pop = 0
    for idx in range(len(d)):
        if d[idx].position <= step:
            count_ele_to_pop += 1
        else:
            break
    
    #pop from left
    for i in range(count_ele_to_pop):
        d.popleft()


def _construct_matrix(d):
    """
    Helper function takes in the window deque and constructs the haplo matrix

    d (deque of RecordLoader objects): window deque

    Returns list of lists: each sublist is a site and each sublist element is the haplotype at that site
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
    Helper counts unique haplotypes in a given window.

    haplo_matrix (list of lists): each sublist is a site and each sublist element is the haplotype at that site

    Returns number of unique haplotypes in the haplo_matrix
    """

    #converting matrix to array and transposing
    matrix = np.array(haplo_matrix)
    matrix = np.transpose(matrix)

    #uniquify
    uniq_matrix = np.unique(matrix, axis=0)

    #counting unique haplotypes
    num_uniq_haps = len(uniq_matrix)


    return num_uniq_haps



### CLASSES ###

class RecordLoader():
    def __init__(self, vcfline):
        """
        This class is for parsing a vcf record for it to be added to a deque.
        The class as two attributes:
        position (int): record position
        GTs (list of str): list of individual genotypes e.g. ["0|0", "1|0"...]
        """

        self.position = int(vcfline.strip().split("\t")[1])
        self.GTs = vcfline.strip().split("\t")[9:]


class ArgLoader():
    def __init__(self, vcfgz_file, windowsize, windowstep):
        """
        Class used to store arguments for the windowed_UPGMA_scan function which is run using the run_hapcount_scan main function.
        This is done so that the main function takes a single argument so that it can be run in parallel.

        vcfgz_file (str): filename or path to file for the vcf to scan
        windowsize (int): SNP window size
        windowstep (int): SNP window step, recommended to be 10% of windowsize
        """

        self.vcf = vcfgz_file
        self.winsize = windowsize
        self.winstep = windowstep

    def __str__(self):
        return self.vcf + "\t" + str(self.winsize) + "\t" + str(self.winstep)




### SLIDING WINDOW SCAN FUNCTION ###

def BPwindow_hap_counter(vcfgz, window_size, window_step):
    """
    Function counts unique haplotypes in sliding windows on a gizped vcf. See other functions for details.
    
    vcfgz (str): gziped vcf file name or path to file.
    SNPwindow_size (int): number of SNPs defining the window size.
    SNPwindow_step (int): number of SNPs to slide window, recommended to be 10% of SNPwindow_size.
    """

    #setting up deque for sliding window
    window_deque = deque()

    #open files
    chrom = vcfgz.split("_")[-1].replace(".vcf.gz", "")
    with gzip.open(vcfgz, "rt") as vcf, open(f"{chrom}_BPwindow{window_size}_BPstep{window_step}_hap_counts.csv", "a") as outcsv:

        #write header line
        outcsv.write(",".join(["CHROM", "START", "END", "NUM_uniq_haps", "SNP_density"]) + "\n")

        #window book keeping
        win_start = 1
        win_end = window_size
        win_step = window_step

        #loop through lines
        for line in vcf:
            #skip header lines
            if line.startswith("#"):
                continue

            ###sliding window scan###

            #holding on  current record
            current_record = RecordLoader(line)

            if current_record.position < win_start:
                 exit("Record position before window start position")
            
            if current_record.position <= win_end:
                window_deque.append(current_record)
                continue


            #increment window until position is reached
            for i in range(win_end, current_record.position, window_step):

                if len(window_deque) == 0:
                    #write to file
                    outcsv.write(_joinany(",", [chrom, win_start, win_end, 0, 0]) + "\n")
                else:
                    #construct haplo matrix
                    hap_m = _construct_matrix(window_deque)
                    #count haps
                    num_haps = _count_unique_haps(hap_m)
                    #compute SNP density
                    snpden = len(window_deque) / window_size
                    #sliding to next step
                    _pop_bp_step(window_deque, win_step)
                    #write to file
                    outcsv.write(_joinany(",", [chrom, win_start, win_end, num_haps, snpden]) + "\n")

                #increment window
                win_end = i + window_step
                win_start = win_end - window_size + 1
                win_step = win_start + window_step

            #append current record to deque
            window_deque.append(current_record)

        #clean up last window
        if len(window_deque) > 0:
            
            #construct haplo matrix
            hap_m = _construct_matrix(window_deque)
            #count haps
            num_haps = _count_unique_haps(hap_m)
            #find number of SNPs
            snpden = len(window_deque) / (win_end - win_start)
            #write to file
            outcsv.write(_joinany(",", [chrom, win_start, win_end, num_haps, snpden]) + "\n")




### MAIN FUNCTION ###
def run_hapcount_scan(argloader_obj):
    """Main Function that runs windowed_UPGMA_scan with a single argument of the class ArgLoader"""

    BPwindow_hap_counter(argloader_obj.vcf, argloader_obj.winsize, argloader_obj.winstep)



# args = ArgLoader("test1.vcf.gz", 100, 10)
# run_hapcount_scan(args)


### PARALLELIZATION ###

def main():

    vcflist = [
        "beagle_phased_biallelic_SNPs_1000GP30X_chr2.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr3.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr4.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr5.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr6.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr7.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr11.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr13.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr14.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr16.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_PARchrX.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_nonPARchrX.vcf.gz", "beagle_imputed_males_biallelic_SNPs_1000GP30X_chrY.vcf.gz"
        ]
    
    loaderlist = []

    for vcffile in vcflist:
        loaderlist.append(ArgLoader(vcffile, 3000, 1500))

    pool = Pool(processes=10)
    pool.map(run_hapcount_scan, loaderlist)


if __name__ == '__main__':
    main()
