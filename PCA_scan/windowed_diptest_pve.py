import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import diptest

#setting output directory
output_dir = os.getcwd()



def diptest_pve_table(filelist, pvefile, dir):
    """
    Loops through a text file of line separated filenames of PCA csv outputs (in that order) and
    performs the dip test on PC1. The function also used to per-window percent variance explained
    file (pvefile) and generates a table with columns: Window, POS, dip_stat, dip_pval, PC1_pve

    Arguments:
    filelist (str): filename of PCA csv files
    pvefile (str): filename of per window PVEs
    dir (str): directory where the data files are contained. This file should be in the current directory!

    Returns DataFrame
    """

    #grabbing list of input PCA csv files
    with open(filelist, "r") as files:
        PCAcsvs = files.readlines()
        PCAcsvs = [x.strip() for x in PCAcsvs]

    #changing working directory
    os.chdir(dir)

    #opening PVE file
    pve_df = pd.read_csv(pvefile)

    #setting lists for new table output
    window = pve_df["Window"].to_list()
    position = pve_df["POS"].to_list()
    dip_stat = []
    dip_pval = []
    pc1pve = pve_df["PC1_variance"].to_list()

    #checking files and pve csv are aligned
    assert len(window) == len(PCAcsvs), "You most likely have selected the wrong files. Windows in pvefile do not match the PCA csvs."

    #looping through windowed PCAs
    for idx, csv in enumerate(PCAcsvs):

        #check for correct window
        win_num = int(csv.split("_")[0].replace("win", ""))
        assert win_num == window[idx], "Windows do not align between pvefile and window PCA."

        #open PCA csv
        winPCA = pd.read_csv(csv)
        pc1 = winPCA["PC1"].to_list()

        #diptest
        dipstat, dipP = diptest.diptest(np.array(pc1))

        #appending to list
        dip_stat.append(dipstat)
        dip_pval.append(dipP)

    #create output df
    df = pd.DataFrame({"Window":window, "POS":position, "Dip_stat":dip_stat, "Dip_pval":dip_pval, "PC1_PVE":pc1pve})


    return df



def plot_dip_pve(dip_pve_df, inversion, title):
    """
    Plots the dip statistic and PVE for PC1 from the DataFrame returned from diptest_pve_table()
    If there is a known inversion in this region, it can be highlighted by setting the inversion argument to
    a tuple of 2 breakpoints. If there is no inversion set to None.
    """

    #creating plot
    fig, axs = plt.subplots(2, 1, figsize=(15, 5))

    #plot dip stat
    axs[0].plot(dip_pve_df["POS"], dip_pve_df["Dip_stat"], color='blue')
    if inversion:
        axs[0].axvspan(inversion[0], inversion[1], color='yellow', alpha=0.4)
    axs[0].set_title("Dip-Test Statistic")
    axs[0].set_ylabel("Dip")
    axs[0].ticklabel_format(style='plain')
    axs[0].set_ylim(0,0.25)

    # #plot dip p-value
    # axs[1].plot(dip_pve_df["POS"], dip_pve_df["Dip_pval"], color='green')
    # axs[1].axhline(y=0.05, color='gray', linestyle='--', linewidth=1)
    # if inversion:
    #     axs[1].axvspan(inversion[0], inversion[1], color='yellow', alpha=0.4)
    # axs[1].set_title("Dip-Test P-value")
    # axs[1].set_ylabel("p-value")

    #plot pve
    axs[1].plot(dip_pve_df["POS"], dip_pve_df["PC1_PVE"], color='red')
    if inversion:
        axs[1].axvspan(inversion[0], inversion[1], color='yellow', alpha=0.4)
    axs[1].set_title("PC1 Percent Variance Explained")
    axs[1].set_ylabel("pve")
    axs[1].set_xlabel("Position")
    axs[1].ticklabel_format(style='plain')
    axs[1].set_ylim(0,1)

    plt.suptitle(title)

    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{title}.png")
    plt.clf()




def main_func(files, pvecsv, directory, plottitle, inv=None):

    newdf = diptest_pve_table(files, pvecsv, directory)
    os.chdir(output_dir)
    plot_dip_pve(newdf, inv, plottitle)





#opening lists of lists
with open("pca_csv_list_of_lists.txt", "r") as lsofls:
    listlist = lsofls.readlines()
    listlist = [x.strip() for x in listlist]

with open("pve_files_list.txt", "r") as pvels:
    pvelist = pvels.readlines()
    pvelist = [x.strip() for x in pvelist]

assert len(listlist) == len(pvelist), "Inputs do not match."


inv_map = {
    "chr9": (113105263, 113119310),
    "chr19": (21647331, 22062459),
    "chr17": (36357258, 37957978),
    "chr15.1": (30077909, 32607507),
    "chr15.2": (22770522, 28852548)
}

for idx in range(len(listlist)):

    chrom = listlist[idx].split("_")[1]
    window = listlist[idx].split("_")[0].replace("win", "")
    inverted = listlist[idx].split("_")[2]
    if inverted == "inv":
        inv_coords = inv_map[chrom]
    else:
        inv_coords = None



    main_func(
        f"pca_csv_file_lists/{listlist[idx]}",
        pvelist[idx],
        f"/Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/sliding_PCA/pca_output/maf0.01_filtered/{window}SNP_windows/{chrom}/{inverted}/",
        f"{chrom} {inverted} in {window} SNP windows",
        inv=inv_coords)
    
