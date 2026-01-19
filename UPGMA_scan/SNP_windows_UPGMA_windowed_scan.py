import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from collections import deque
from multiprocessing import Pool
import sys
sys.setrecursionlimit(10000)


### HELPER FUNCTIONS ###

def _joinany(sep, mylist):
    """Function to join a list by a sep regardless of element types"""

    strs = [str(x) for x in mylist]

    return sep.join(strs)


# def _plot_dendro(linkage_array, title):
#     """
#     Helper function used to plot selected tree and save to png. 
#     It is only called within specific intervals (see windowed_UPGMA_scan)
#     """

#     dendrogram(linkage_array)
#     plt.xticks([])
#     plt.title(title.replace("_", " "))
#     plt.savefig(f"{title}.png")




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

    d (deque of RecordLoader objects): window deque

    Returns list of lists: each sublist is a site and each sublist element is the haplotype at that site
    """

    matrix = []
    for record in d:
        hap_ls = []
        for gt in record.GTs:
            # if "|" in gt:
            #     haps = gt.split("|")
            #     hap_ls.append(int(haps[0]))
            #     hap_ls.append(int(haps[1]))
            # else:
            #     hap_ls.append(int(gt))
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

    return avg_branch



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
    
    #finding average branch length
    avglen = _find_avg_branch(tre)


    return avglen       #, tre



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
    def __init__(self, vcfgz_file, windowsize, windowstep):     #, intervals_to_plot=None
        """
        Class used to store arguments for the windowed_UPGMA_scan function which is run using the run_UPGMA_scan main function.
        This is done so that the main function takes a single argument so that it can be run in parallel.

        vcfgz_file (str): filename or path to file for the vcf to scan
        windowsize (int): SNP window size
        windowstep (int): SNP window step, recommended to be 10% of windowsize
        #####intervals_to_plot (default is None, set to [(int,int,int),...]): genomic intervals (1st and 2nd int) in which to plot some trees. The 3rd int specifies to plot a tree every n windows.
        """

        self.vcf = vcfgz_file
        self.winsize = windowsize
        self.winstep = windowstep
        # self.plot_here = intervals_to_plot

    def __str__(self):
        return self.vcf + "\t" + str(self.winsize) + "\t" + str(self.winstep)




### SLIDING WINDOW SCAN FUNCTION ###

def windowed_UPGMA_scan(vcfgz, SNPwindow_size, SNPwindow_step):     #, plotting_intervals
    """
    Function performs the windowed UPGMA scan on a gizped vcf, and outputs 
    Average Branch Length results to a csv. See other functions for details.
    
    vcfgz (str): gziped vcf file name or path to file.
    SNPwindow_size (int): number of SNPs defining the window size.
    SNPwindow_step (int): number of SNPs to slide window, recommended to be 10% of SNPwindow_size.
    #####plotting_intervals ([(int, int, int)] or None): genomic position intervals (1st and 2nd int) in which to plot UPGMA trees. The 3rd int specifies to plot a tree every n windows.
    """

    #setting up deque for sliding window
    window_deque = deque()

    #window counter
    counter = -1

    #open files
    chrom = vcfgz.split("_")[-1].replace(".vcf.gz", "")
    with gzip.open(vcfgz, "rt") as vcf, open(f"{chrom}_window{SNPwindow_size}_step{SNPwindow_step}_avg_branch_len.csv", "a") as outcsv:

        #write header line
        outcsv.write(",".join(["CHROM", "START", "END", "AVG_branch_length"]) + "\n")

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

                    #progress window counter (0-based index)
                    counter += 1

                    #grab window position
                    start_pos = window_deque[0].position
                    end_pos = window_deque[-1].position + 1
                    # win_pos = (start_pos + end_pos) // 2

                    #construct haplo matrix
                    hap_m = _construct_matrix(window_deque)

                    #run UPGMA
                    avgbranchlen = _runUPGMA(hap_m)     #, tree_array

                    # #plot tree?
                    # if plotting_intervals:
                    #     for start, end, mod_num in plotting_intervals:
                    #         if win_pos >= start and win_pos <= end and (counter % mod_num) == 0:
                    #             tree_title = f"window{counter}_{chrom}_from_{start_pos}_to_{end_pos}"
                    #             _plot_dendro(tree_array, tree_title)
                    #             break


                    #sliding to next step
                    _pop_step(window_deque, SNPwindow_step)

                    #write to file
                    outcsv.write(_joinany(",", [chrom, start_pos, end_pos, avgbranchlen]) + "\n")

                else:
                    continue

        #clean up last window
        if len(window_deque) > 0:
            
            #grab window position
            start_pos = window_deque[0].position
            end_pos = window_deque[-1].position + 1
            # win_pos = (window_deque[0].position + window_deque[-1].position) // 2

            #construct haplo matrix
            hap_m = _construct_matrix(window_deque)

            #run UPGMA
            avgbranchlen = _runUPGMA(hap_m)     #, tree_array
            
            #write to file
            outcsv.write(_joinany(",", [chrom, start_pos, end_pos, avgbranchlen]) + "\n")




### MAIN FUNCTION ###
def run_UPGMA_scan(argloader_obj):
    """Main Function that runs windowed_UPGMA_scan with a single argument of the class ArgLoader"""

    windowed_UPGMA_scan(argloader_obj.vcf, argloader_obj.winsize, argloader_obj.winstep)    #, argloader_obj.plot_here



# args = ArgLoader("test1.vcf.gz", 100, 10)
# run_UPGMA_scan(args)



### PARALLELIZATION ###

def main():
    if len(sys.argv) < 2:
        exit("missing argument")
    threads = int(sys.argv[1])

    # plot_in_these_intervals = {
    #     "chr3":[(195622595, 197667517, 250), (232674321, 233674321, 250)],
    #     "chr5":[(64466001, 64482016, 5), (43268583, 43278583, 5)],
    #     "chr8":[(7064966, 12716088, 250), (42683243, 43683243, 250)],
    #     "chr11":[(71564791, 71583752, 5), (93689963, 93699963, 5)],
    #     "chr15":[(30077909, 32607507, 250), (33264532, 34264532, 250)],
    #     "chr16":[(75205642, 75223319, 5), (53674896, 53684896, 5)],
    #     "chr17":[(45495836, 46707123, 250), (48435128, 49435128, 250)]
    # }
    vcflist = [
        "beagle_phased_biallelic_SNPs_1000GP30X_chr2.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr3.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr4.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr5.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr6.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr7.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr11.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr13.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr14.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_chr16.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz",
        "beagle_phased_biallelic_SNPs_1000GP30X_PARchrX.vcf.gz", "beagle_phased_biallelic_SNPs_1000GP30X_nonPARchrX.vcf.gz", "beagle_imputed_males_biallelic_SNPs_1000GP30X_chrY.vcf.gz"
        ]
    windowing = [(1000, 500), (10000, 5000), (100, 50)]    #(100, 10), 
    loaderlist = []

    for size, step in windowing:
        for vcffile in vcflist:
            # chromosome = vcffile.split("_")[-1].replace(".vcf.gz", "")
            # if chromosome in plot_in_these_intervals:
            #     loaderlist.append(ArgLoader(vcffile, size, step, intervals_to_plot=plot_in_these_intervals[chromosome]))
            # else:
            loaderlist.append(ArgLoader(vcffile, size, step))

    pool = Pool(processes=threads)
    pool.map(run_UPGMA_scan, loaderlist)


if __name__ == '__main__':
    main()
