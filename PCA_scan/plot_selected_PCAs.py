import pandas as pd
import matplotlib.pyplot as plt



def plot6scatters(listfile, inversion_samples):
    """
    Function takes a file with a line separated list of 6 paths to PCA csvs
    and plots the scatter plots in a 3x2 grid.
    """
    #grab main title string
    title = listfile.replace(".txt", "").replace("_", " ")

    #check in inversion or control
    if "inv" in title:
        inverted = True
    else:
        inverted = False

    #open list of paths file
    with open(listfile, "r") as listofpaths:
        paths = listofpaths.readlines()
        paths = [x.strip() for x in paths]

    #setup plot
    fig, axs = plt.subplots(2, 3, figsize=(10, 8))
    if inverted:
        point_color = "red"
    else:
        point_color = "blue"
    plt.suptitle(title)

    #loop through csvs
    for idx, csv in enumerate(paths):
        
        #find subplot coordinates
        if idx >= 3:
            row = 1
            col = idx - 3
        else:
            row = 0
            col = idx

        #grab subtitle
        subtitle = csv.split("/")[-1].split("_")[0].replace("win", "Window ")

        #open table
        df = pd.read_csv(csv)

        #plotting PCs
        axs[row, col].scatter(df["PC1"], df["PC2"], color=point_color, s=1)
        if inversion_samples:
            for idx in inversion_samples:
                axs[row, col].scatter(df["PC1"][idx], df["PC2"][idx], color="y", s=20, alpha=0.8)
        axs[row, col].set_title(subtitle)
        axs[row, col].set_xlabel("PC1")
        axs[row, col].set_ylabel("PC2")

    plt.tight_layout()
    #plt.show()
    plt.savefig(f"{title}.png")
    plt.clf()






listoflistfiles = [
    "25SNPwin_chr19_inv.txt",
    "25SNPwin_chr17_inv.txt",
    "25SNPwin_chr15.1_inv.txt",
    "25SNPwin_chr15.2_inv.txt",
    "50SNPwin_chr19_inv.txt",
    "50SNPwin_chr17_inv.txt",
    "50SNPwin_chr15.1_inv.txt",
    "50SNPwin_chr15.2_inv.txt",
    "25SNPwin_chr19_control.txt",
    "25SNPwin_chr17_control.txt",
    "25SNPwin_chr15.1_control.txt",
    "25SNPwin_chr15.2_control.txt",
    "50SNPwin_chr19_control.txt",
    "50SNPwin_chr17_control.txt",
    "50SNPwin_chr15.1_control.txt",
    "50SNPwin_chr15.2_control.txt"
]



#open sample files
with open("1000GP_samples.txt", "r") as samplesf:
    samples = samplesf.readlines()
    samples = [x.strip() for x in samples]

with open("chr19_inv_samples.txt", "r") as chr19samplesf:
    chr19_inv_samples = chr19samplesf.readlines()
    chr19_inv_samples = [x.strip() for x in chr19_inv_samples]

with open("chr9_inv_samples.txt", "r") as chr9samplesf:
    chr9_inv_samples = chr9samplesf.readlines()
    chr9_inv_samples = [x.strip() for x in chr9_inv_samples]

with open("chr17_inv_samples.txt", "r") as chr17samplesf:
    chr17_inv_samples = chr17samplesf.readlines()
    chr17_inv_samples = [x.strip() for x in chr17_inv_samples]

with open("chr15.1_inv_samples.txt", "r") as chr15samplesf:
    chr15_inv_samples = chr15samplesf.readlines()
    chr15_inv_samples = [x.strip() for x in chr15_inv_samples]


def find_inv_sample_idxs(all_samples, inv_samples):
    """
    Function to find which idxs the inversion samples are at in the all_samples list.
    Note that by "idxs" I mean the idx of the haplotypes coresponding to that sample.
    e.g. sample idx of 2 is haplotype idxs of 4 and 5
    """
    haplo_idxs = []
    for idx, sample in enumerate(all_samples):
        for inv_s in inv_samples:
            if inv_s == sample:
                haplo_idxs.append(idx*2)
                haplo_idxs.append((idx*2)+1)
                break
            else:
                continue

    return haplo_idxs


inv_sample_idx_map = {}
inv_sample_idx_map["chr19"] = find_inv_sample_idxs(samples, chr19_inv_samples)
inv_sample_idx_map["chr9"] = find_inv_sample_idxs(samples, chr9_inv_samples)
inv_sample_idx_map["chr15.1"] = find_inv_sample_idxs(samples, chr15_inv_samples)
inv_sample_idx_map["chr17"] = find_inv_sample_idxs(samples, chr17_inv_samples)



for lsfile in listoflistfiles:

    chrom = lsfile.split("_")[1]
    if chrom in inv_sample_idx_map.keys():
        inv_s_idxs = inv_sample_idx_map[chrom]
    else:
        inv_s_idxs = None

    plot6scatters(lsfile, inv_s_idxs)
    