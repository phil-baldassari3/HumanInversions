import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from multiprocessing import Pool




def _get_sitehaplist(record):
    """
    Helper function that take in a line (record) from vcf and Returns a site haplotype list needed for appending to matrix
    """

    GT_ls = record.strip().split("\t")[9:]
    hap_ls = []
    for gt in GT_ls:
        haps = gt.split("|")
        hap_ls.append(int(haps[0]))
        hap_ls.append(int(haps[1]))

    return hap_ls


def _get_record_position(record):
    """
    Helper function grabs the SNP coordinate for the input record and returns it as an integer
    """

    record_pos = int(record.strip().split("\t")[1])

    return record_pos




def find_top_branch(linkage_array):
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

    Returns top_branch, norm_top_branch
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


    return top_branch, norm_top_branch




# def plot6dendros(list_of_arrays, plottitle):
#     """
#     Function takes a list of linakge arrays and plots the dendrograms in a 3x2 grid.
#     """


#     #setup plot
#     fig, axs = plt.subplots(2, 3, figsize=(10, 8))
#     plt.suptitle(plottitle)

#     #loop through csvs
#     for idx, tree in enumerate(list_of_arrays):
        
#         #find subplot coordinates
#         if idx >= 3:
#             row = 1
#             col = idx - 3
#         else:
#             row = 0
#             col = idx

#         #plotting
#         dendrogram(tree, ax=axs[row, col])

#     plt.tight_layout()
#     #plt.show()
#     plt.savefig(f"{plottitle}.png")
#     plt.clf()


    


def UPGMA_on_SNP_windows(vcfgz, SNPwindow, output_dir):
    """
    Runs UPGMA hierarchical clustering on haplotypes from a phased vcf in non-overlapping SNP windows

    The function will slide a window of SNPwindow size through the vcf, extract haplotypes, and run
    UPGMA using scipy.cluster.hierarch.linkage. For each window, a text file containing the linkage array
    output will be saved to the output_dir. It will then output a csv file of top branch length for 
    each window will be saved to the current dir. It also saves a linakge array as a .npy for each window
    to the output_dir
    """

    #empty lists for variance explained df
    win = []
    pos = []
    longest_branch = []
    norm_longest_branch = []

    #open file
    with gzip.open(vcfgz, "rt") as vcf:

        #setting counters
        window_counter = 0

        SNP_counter = 0
        start_pos = 0
        end_pos = 0

        #empty matrix
        matrix = []

        #loop through lines
        for line in vcf:

            #skip header
            if line.startswith("#"):
                continue

            #run through SNPs
            else:
                
                #update counters
                if SNP_counter == 0:
                    start_pos = _get_record_position(line)
                SNP_counter += 1
                
                #adding record to hap matrix
                site_hap_ls = _get_sitehaplist(line)
                matrix.append(site_hap_ls)

                #if end of window is reached
                if SNP_counter == SNPwindow:
                    #logging meta data
                    win.append(window_counter)
                    end_pos = _get_record_position(line)
                    pos.append((start_pos + end_pos) // 2)

                    #converting matrix to array and transposing
                    matrix = np.array(matrix)
                    matrix = np.transpose(matrix)

                    #running UPGMA
                    tre = linkage(matrix, method="average", metric="euclidean")
                    longest, norm_longest = find_top_branch(tre)
                    longest_branch.append(longest)
                    norm_longest_branch.append(norm_longest)
                    np.save(f"{output_dir}win{window_counter}_windowsize{SNPwindow}_from_{start_pos}_to_{end_pos}_{vcfgz}".replace(".vcf.gz", ".npy"), tre)


                    #resetting counts
                    SNP_counter = 0
                    start_pos = 0
                    end_pos = 0

                    #counting window
                    window_counter += 1

                    #clear matrix
                    matrix = []

                else:
                    continue

        #cleanup
        if SNP_counter != 0:
            win.append(window_counter)
            end_pos = _get_record_position(line)
            pos.append((start_pos + end_pos) // 2)

            matrix = np.array(matrix)
            matrix = np.transpose(matrix)

            tre = linkage(matrix, method="average", metric="euclidean")
            longest, norm_longest = find_top_branch(tre)
            longest_branch.append(longest)
            norm_longest_branch.append(norm_longest)
            np.save(f"{output_dir}win{window_counter}_windowsize{SNPwindow}_from_{start_pos}_to_{end_pos}_{vcfgz}".replace(".vcf.gz", ".npy"), tre)

            #clear matrix
            matrix = []

    #saving variance explained csv
    vdf = pd.DataFrame({"Window":win, "POS":pos, "Longest_branch_length":longest_branch, "Normalized_longest_branch_length":norm_longest_branch})
    vdf.to_csv(f"longest_branch_lengths{SNPwindow}_{vcfgz}".replace(".vcf.gz", ".csv"), index=False)

    print(f"Done sliding window UPGMA for {vcfgz} with window size {SNPwindow}")





class ClustFuncLoader():
    def __init__(self, vcfgzfile, window):
        """
        Class used to load arguments into the main_func so that it can be used with multiprocessing
        """

        self.filename = vcfgzfile
        self.window_size = window

        chrom = vcfgzfile.split("_")[0]

        if "inv" in vcfgzfile:
            inverted = "inv"
        else:
            inverted = vcfgzfile.split("_")[1]

        self.output_directory = f"clustering_output/{window}SNP_windows/{chrom}/{inverted}/"




def main_func(clustfuncloader_obj):
    """
    Function runs UPGMA_on_SNP_windows() but takes in a single argument which is of the ClustFuncLoader type
    """
    UPGMA_on_SNP_windows(clustfuncloader_obj.filename, clustfuncloader_obj.window_size, clustfuncloader_obj.output_directory)




#files to run on
# vcflist = [
#     "chr15.1_1MB_inv_1MB_4529599bp.vcf.gz", "chr15.2_1MB_inv_1MB_8082027bp.vcf.gz", "chr15.1_control1_4529599bp.vcf.gz",
#     "chr15.2_control1_8082027bp.vcf.gz", "chr15.1_control2_4529599bp.vcf.gz", "chr15.2_control2_8082027bp.vcf.gz",
#     "chr17_1MB_inv_1MB_3600721bp.vcf.gz", "chr17_control1_3600721bp.vcf.gz", "chr17_control2_3600721bp.vcf.gz",
#     "chr19_1MB_inv_1MB_2415129bp.vcf.gz", "chr19_control1_2415129bp.vcf.gz", "chr19_control2_2415129bp.vcf.gz",
#     "chr9_1MB_inv_1MB_2014048bp.vcf.gz", "chr9_control1_2014048bp.vcf.gz", "chr9_control2_2014048bp.vcf.gz"
# ]

vcflist = [
    "chr15.1_1MB_inv_1MB_4529599bp.vcf.gz", "chr15.2_1MB_inv_1MB_8082027bp.vcf.gz", "chr15.1_control1_4529599bp.vcf.gz",
    "chr15.2_control1_8082027bp.vcf.gz", "chr15.1_control2_4529599bp.vcf.gz", "chr15.2_control2_8082027bp.vcf.gz",
    "chr17_1MB_inv_1MB_3600721bp.vcf.gz", "chr17_control1_3600721bp.vcf.gz", "chr17_control2_3600721bp.vcf.gz"
]


#initializing ClustFuncLoader objects
loaderlist = []
for win in [200, 400]:
    for file in vcflist:
        loader = ClustFuncLoader(file, win)
        loaderlist.append(loader)





#running main_func in parallel with arguments from ClustFuncLoader
def run_in_parallel():
    pool = Pool(processes=5)
    pool.map(main_func, loaderlist)


if __name__ == '__main__':
    run_in_parallel()