import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import deque
from multiprocessing import Pool
import sys


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



def _find_avg_branch(linkage_array):
    """Helper function that finds the noramlized average branch length from a linakge array"""

    #tree height
    tree_height = linkage_array[-1, 2]

    #lengths
    lengths = []
    for clust in linkage_array:
        lengths.append(clust[2])

    #averaging
    avg_branch = (sum(lengths) / len(lengths)) / tree_height

    return avg_branch, tree_height


def _find_top_branch(linkage_array):
    """
    The scipy.cluster.hierarchy,linkage array is confusing. Luckily, Dr. Jim knew exactly how to parse it
    This function finds the longest branch length (which is not the same as the maximum value in the branch
    column of the linakge output array) by first taking the maximum length value in the array x 2 (this is 
    the height of the tree). Next, this row of the array contains the indeces representing the 2 child clusters
    except these indeces may be greater than the length of the number of data points if they are superclusters
    of smaller clusters. Therefore in order to subtract the branch length of these child clusters from the tree
    height, we find the indeces for these child clusters by taking the supercluster index - number of original
    data points. We subtract these children branch lengths from the tree height to find the top branch length. 
    The branch length is also then normalized by the height of the tree.

    e.g.
    data = [[0 0 0 0 0 0 0 0 0 0]
            [1 1 1 1 1 1 1 1 1 1]
            [1 1 1 0 1 1 1 1 1 1]
            [1 1 1 1 1 1 1 1 1 1]]

    linkage_array = [   [1.         3.         0.         2.        ]
                        [2.         4.         1.         3.        ]
                        [0.         5.         3.10818511 4.        ]   ]

    longest_branch = (2*3.10818511) - 
                        (linkage_array[(linkage_array[-1,0]-linkage_array[-1, 3],2]) + 
                        (linkage_array[(linkage_array[-1,1]-linkage_array[-1, 3],2]))
                        if linkage_array[-1,0]-linkage_array[-1, 3] >= 0, else subtract 0
                    = 6.216 - (0 + 1) = 5.216
    norm_longest_branch = longest_branch / 2*3.10818511 = 0.839

    Returns norm_top_branch
    """

    #finding tree height
    tree_height = linkage_array[-1, 2] * 2

    #finding number of data points (leaves)
    n = linkage_array[-1, 3]

    #finding children heights
    child1idx = int(linkage_array[-1, 0] - n)
    child2idx = int(linkage_array[-1, 1] - n)

    if child1idx >= 0:
        child1height = linkage_array[child1idx, 2]
    else:
        child1height = 0

    if child2idx >= 0:
        child2height = linkage_array[child2idx, 2]
    else:
        child2height = 0

    #computing longest branch length
    top_branch = tree_height - (child1height + child2height)
    norm_top_branch = top_branch / tree_height


    return norm_top_branch



def _runUPGMA(haplo_matrix):
    """
    Helper function runs UPGMA on haplotypes.

    haplo_matrix (list of lists): each sublist is a site and each sublist element is the haplotype at that site

    Returns float for normalized average branch length #####and ndarray from linkage()
    """

    #converting matrix to array and transposing
    matrix = np.array(haplo_matrix)
    matrix = np.transpose(matrix)

    #running UPGMA
    tre = linkage(matrix, method="average", metric="euclidean")
    
    #finding average branch length, tree height, and longest branch
    avglen, tree_height = _find_avg_branch(tre)
    longestbranch = _find_top_branch(tre)


    return avglen, longestbranch, tree_height



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
        Class used to store arguments for the windowed_UPGMA_scan function which is run using the run_UPGMA_scan main function.
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

def windowed_UPGMA_scan(vcfgz, window_size, window_step):
    """
    Function performs the windowed UPGMA scan on a gizped vcf, and outputs 
    Average Branch Length results to a csv. See other functions for details.
    
    vcfgz (str): gziped vcf file name or path to file.
    SNPwindow_size (int): number of SNPs defining the window size.
    SNPwindow_step (int): number of SNPs to slide window, recommended to be 10% of SNPwindow_size.
    """

    #setting up deque for sliding window
    window_deque = deque()

    #open files
    chrom = vcfgz.split("_")[-1].replace(".vcf.gz", "")
    with gzip.open(vcfgz, "rt") as vcf, open(f"{chrom}_window{window_size}_step{window_step}_avg_branch_len.csv", "a") as outcsv:

        #write header line
        outcsv.write(",".join(["CHROM", "START", "END", "AVG_branch_length", "LONGEST_branch_length", "Tree_Height", "SNP_density"]) + "\n")

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
                    outcsv.write(_joinany(",", [chrom, win_start, win_end, 0, 0, 0, 0]) + "\n")
                else:
                    #construct haplo matrix
                    hap_m = _construct_matrix(window_deque)
                    #run UPGMA
                    avgbranchlen, longest_branch_len, height_of_tree = _runUPGMA(hap_m)
                    #find SNP density
                    snpden = len(window_deque) / window_size
                    #sliding to next step
                    _pop_bp_step(window_deque, win_step)
                    #write to file
                    outcsv.write(_joinany(",", [chrom, win_start, win_end, avgbranchlen, longest_branch_len, height_of_tree, snpden]) + "\n")

                #increment window
                win_end = i + window_step
                win_start = win_end - window_size + 1
                win_step = win_start + window_step

                # win_start += window_step
                # win_end += window_step
                # win_step += window_step

            #append current record to deque
            window_deque.append(current_record)



        #clean up last window
        if len(window_deque) > 0:
            
            #construct haplo matrix
            hap_m = _construct_matrix(window_deque)
            #run UPGMA
            avgbranchlen, longest_branch_len, height_of_tree = _runUPGMA(hap_m)
            #find SNP density
            snpden = len(window_deque) / (win_end - win_start)
            #write to file
            outcsv.write(_joinany(",", [chrom, win_start, win_end, avgbranchlen, longest_branch_len, height_of_tree, snpden]) + "\n")




### MAIN FUNCTION ###
def run_UPGMA_scan(argloader_obj):
    """Main Function that runs windowed_UPGMA_scan with a single argument of the class ArgLoader"""

    windowed_UPGMA_scan(argloader_obj.vcf, argloader_obj.winsize, argloader_obj.winstep)



# args = ArgLoader("beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz", 3000, 1500)
# run_UPGMA_scan(args)


### PARALLELIZATION ###

def main():

    # vcflist = [
    #     "beagle_phased_biallelic_SNPs_1000GP30X_chr2.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr3.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr4.vcf.gz",
    #     "beagle_phased_biallelic_SNPs_1000GP30X_chr5.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr6.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr7.vcf.gz",
    #     "beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr11.vcf.gz",
    #     "beagle_phased_biallelic_SNPs_1000GP30X_chr13.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr14.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz",
    #     "beagle_phased_biallelic_SNPs_1000GP30X_chr16.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz",
    #     "beagle_phased_biallelic_SNPs_1000GP30X_PARchrX.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_nonPARchrX.vcf.gz", "beagle_imputed_males_biallelic_SNPs_1000GP30X_chrY.vcf.gz"
    #     ]

    vcflist = [
        "beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz", 
        "beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz",
        ]
    
    loaderlist = []

    for vcffile in vcflist:
        loaderlist.append(ArgLoader(vcffile, 3000, 1500))

    pool = Pool(processes=4)
    pool.map(run_UPGMA_scan, loaderlist)


if __name__ == '__main__':
    main()
