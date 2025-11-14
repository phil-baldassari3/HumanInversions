import os
import glob
import pandas as pd
import numpy as np


def _branch_length(linkage_array, cluster_idx, n):
    """Helper function to find branch length of a cluster (row) in a linakge array"""

    #find heights
    tree_height = linkage_array[-1, 2] * 2
    cluster_height = linkage_array[cluster_idx, 2] * 2

    #finding children heights
    child1idx = int(linkage_array[cluster_idx, 0] - n)
    child2idx = int(linkage_array[cluster_idx, 1] - n)

    if child1idx >= 0:
        child1height = linkage_array[child1idx, 2]
    else:
        child1height = 0

    if child2idx >= 0:
        child2height = linkage_array[child2idx, 2]
    else:
        child2height = 0

    #computing branch length
    raw_branch_len = cluster_height - (child1height + child2height)
    branch_len = raw_branch_len / tree_height

    return branch_len



def find_longest_branch(linkage_array):
    """
    Function iterates through linkage array in reverse order and finds the branch length
    of each cluster using _branch_length(). After finding all the branch length, this function
    then returns the maximum value.
    """
    #empty branch length list
    branch_lengths = []

    #grab number of tips
    num_tips = linkage_array[-1, 3]

    #iterating through linakge array
    for idx in reversed(range(len(linkage_array))):
        length = _branch_length(linkage_array, idx, num_tips)
        branch_lengths.append(length)

    return max(branch_lengths)
  



def find_avg_branch(linkage_array):
    """
    Finds the noramlized average branch length from a linakge array
    """

    #tree height
    tree_height = linkage_array[-1, 2]

    #lengths
    lengths = []
    for clust in linkage_array:
        lengths.append(clust[2])

    #averaging
    avg_branch = (sum(lengths) / len(lengths)) / tree_height

    return avg_branch



def longest_and_avg_braches(directory):
    """
    Loops through directory containing .npy files and outputs as csv with
    longest branch length and average branch length.
    """

    #how many files in directory
    files = os.listdir(directory)
    files = [x for x in files if ".npy" in x]
    num_files = len(files)
    

    #empty lists for output table
    win = []
    pos = []
    longest_branch = []
    avg_branch = []

    #looping through files in order:
    for i in range(num_files):
        globls = glob.glob(f"{directory}/win{i}_*.npy")
        assert len(globls) == 1, "you globbed wrong!"
        npy = globls[0]
        npy_name = npy.split("/")[-1]
        npy_array = np.load(npy)
        
        #grabbing window and position
        window = int(npy_name.split("_")[0].replace("win", ""))
        win.append(window)

        pos1 = int(npy_name.split("_")[3])
        pos2 = int(npy_name.split("_")[5])
        position = (pos1 + pos2) // 2
        pos.append(position)

        #compute longest branch
        longest = find_longest_branch(npy_array)
        longest_branch.append(longest)

        #compute average branch
        avg = find_avg_branch(npy_array)
        avg_branch.append(avg)

    #making df
    df = pd.DataFrame({"Window":win, "POS":pos, "Longest_branch_length":longest_branch, "Average_branch_length":avg_branch})

    windowing = directory.split("/")[-3]
    chrom = directory.split("/")[-2]
    inverted = directory.split("/")[-1]
    df.to_csv(f"LongestAvg_branch_lengths_{windowing}_{chrom}_{inverted}.csv", index=False)
        





longest_and_avg_braches("clustering_output/200SNP_windows/chr9/inv")
longest_and_avg_braches("clustering_output/200SNP_windows/chr9/control1")
longest_and_avg_braches("clustering_output/200SNP_windows/chr9/control2")

longest_and_avg_braches("clustering_output/200SNP_windows/chr19/inv")
longest_and_avg_braches("clustering_output/200SNP_windows/chr19/control1")
longest_and_avg_braches("clustering_output/200SNP_windows/chr19/control2")

longest_and_avg_braches("clustering_output/200SNP_windows/chr17/inv")
longest_and_avg_braches("clustering_output/200SNP_windows/chr17/control1")
longest_and_avg_braches("clustering_output/200SNP_windows/chr17/control2")

longest_and_avg_braches("clustering_output/200SNP_windows/chr15.1/inv")
longest_and_avg_braches("clustering_output/200SNP_windows/chr15.1/control1")
longest_and_avg_braches("clustering_output/200SNP_windows/chr15.1/control2")

longest_and_avg_braches("clustering_output/200SNP_windows/chr15.2/inv")
longest_and_avg_braches("clustering_output/200SNP_windows/chr15.2/control1")
longest_and_avg_braches("clustering_output/200SNP_windows/chr15.2/control2")



longest_and_avg_braches("clustering_output/400SNP_windows/chr9/inv")
longest_and_avg_braches("clustering_output/400SNP_windows/chr9/control1")
longest_and_avg_braches("clustering_output/400SNP_windows/chr9/control2")

longest_and_avg_braches("clustering_output/400SNP_windows/chr19/inv")
longest_and_avg_braches("clustering_output/400SNP_windows/chr19/control1")
longest_and_avg_braches("clustering_output/400SNP_windows/chr19/control2")

longest_and_avg_braches("clustering_output/400SNP_windows/chr17/inv")
longest_and_avg_braches("clustering_output/400SNP_windows/chr17/control1")
longest_and_avg_braches("clustering_output/400SNP_windows/chr17/control2")

longest_and_avg_braches("clustering_output/400SNP_windows/chr15.1/inv")
longest_and_avg_braches("clustering_output/400SNP_windows/chr15.1/control1")
longest_and_avg_braches("clustering_output/400SNP_windows/chr15.1/control2")

longest_and_avg_braches("clustering_output/400SNP_windows/chr15.2/inv")
longest_and_avg_braches("clustering_output/400SNP_windows/chr15.2/control1")
longest_and_avg_braches("clustering_output/400SNP_windows/chr15.2/control2")

