import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import sys
sys.setrecursionlimit(10000)



def plot6dendros(list_of_npys):
    """
    Function takes a file listing of npy array files and plots the dendrograms in a 3x2 grid.
    """

    #grab main title string
    title = list_of_npys.replace(".txt", "").replace("_", " ")

    #open list of paths file
    with open(list_of_npys, "r") as listofpaths:
        paths = listofpaths.readlines()
        paths = [x.strip() for x in paths]

    #setup plot
    fig, axs = plt.subplots(2, 3, figsize=(15, 5))

    #loop through csvs
    for idx, npy in enumerate(paths):

        #grab subtitle
        subtitle = npy.split("/")[-1].split("_")[0].replace("win", "Window ")
    
        #find subplot coordinates
        if idx >= 3:
            row = 1
            col = idx - 3
        else:
            row = 0
            col = idx

        #plotting
        tree = np.load(npy)
        dendrogram(tree, ax=axs[row, col])
        axs[row, col].set_title(subtitle)
        axs[row, col].set_xticks([])
        #axs[row, col].set_ylim(0,11)

    plt.suptitle(title)
    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{title}.png")
    plt.clf()






listoflistfiles = [
    "200SNPwin_chr19_inv.txt",
    "200SNPwin_chr17_inv.txt",
    "200SNPwin_chr15.1_inv.txt",
    "200SNPwin_chr15.2_inv.txt",
    "400SNPwin_chr19_inv.txt",
    "400SNPwin_chr17_inv.txt",
    "400SNPwin_chr15.1_inv.txt",
    "400SNPwin_chr15.2_inv.txt",
    "200SNPwin_chr19_control.txt",
    "200SNPwin_chr17_control.txt",
    "200SNPwin_chr15.1_control.txt",
    "200SNPwin_chr15.2_control.txt",
    "400SNPwin_chr19_control.txt",
    "400SNPwin_chr17_control.txt",
    "400SNPwin_chr15.1_control.txt",
    "400SNPwin_chr15.2_control.txt"
]





#testing
# p = "clustering_output/200SNP_windows/chr15.1/control1/win431_windowsize200_from_68584965_to_68590110_chr15.1_control1_4529599bp.npy"
# tree = np.load(p)
# dendrogram(tree)
# plt.xticks([])
# plt.show()



for lsfile in listoflistfiles:

    plot6dendros(lsfile)