import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


inv_map = {
    "chr1":[(205209479, 205209706), (187495697, 187497600)],
    "chr4":[(40233407, 40235440), (3894504, 9763113)],
    "chr6":[(130527041, 130531150), (31041443, 31042319)],
    "chr14":[(60604530, 60613248), (34540845, 34562270)],
    "chr15":[(22770522, 28852548), (30077909, 32607507)]
    }


def estimate_heterozygosity_from_af(freqfile, SNPwindow_size, SNPwindow_step):
    """
    Function takes a .freq output from vcftools and estimates
    windowed heterozygosity.

    Returns DataFrame of windowed heterozygosity
    """

    #open file
    cols = ["CHROM", "POS", "N_ALLELES", "N_CHR", "p1", "p2"]
    freqdf = pd.read_csv(freqfile, sep="\t", header=None, skiprows=1, names=cols)

    #clean data
    freqdf["p1"] = freqdf["p1"].apply(lambda x: float(x.split(":")[1]))
    freqdf["p2"] = freqdf["p2"].apply(lambda x: float(x.split(":")[1]))

    #compute per site heterozygosity
    freqdf["Het"] = 1 - ((freqdf["p1"]**2) + (freqdf["p2"]**2))

    #compute windowed heterozygosity
    pos = []
    avgHet = []

    for start in range(0, len(freqdf)-SNPwindow_size, SNPwindow_step):
        stop = start + SNPwindow_size

        window = freqdf[start:stop]
        pos.append(window["POS"].mean())
        avgHet.append(window["Het"].mean())

    windowedHet = pd.DataFrame({"Position":pos, "Heterozygosity":avgHet})

    return windowedHet


def load_pi_data(pifile):
    """
    Function loads in windowed pi data from vcftools.
    THe main use for this function is to take the average b/n
    the window start and stop positions for better plotting

    Returns DataFrame of windowed Pi
    """

    #open file
    pidf = pd.read_csv(pifile, sep="\t")

    #average start and end positions
    pidf["Position"] = (pidf["BIN_START"] + pidf["BIN_END"]) // 2

    #select columns
    windowedpi = pidf[["Position", "PI"]]

    return windowedpi


def plot_stat_w_groundtruth_inv(data, stat, chrom, inversion_map, window):
    """
    Plots either population level windowed Heterozygosity or PI and 
    highlights the spans of groundtruth inversions from a dict map given
    """

    #groundtruth inversions
    grndt_invs = inversion_map[chrom]

    #window labels
    if stat == "Heterozygosity":
        win_lab = "SNP windows"
    elif stat == "PI":
        win_lab = "bp windows"
    else:
        print("oops")

    plt.figure(figsize=(15, 5))
    plt.plot(data["Position"], data[stat])
    plt.xlabel(chrom)
    plt.ylabel(stat)
    plt.ticklabel_format(style='plain')
    
    for span in grndt_invs:
        plt.axvspan(span[0], span[1], color='red', alpha=0.5)

    plt.title(f"{chrom} {stat} {window} {win_lab}")
    #plt.savefig(f"{chrom}_{stat}_{window}_{win_lab}.png")
    plt.show()
    plt.clf()



def main_func(datafile, stat, chromosome, windowsize, windowstep):
    """
    Main function that loads in the data and plots
    """

    #loading data
    if stat == "Heterozygosity":
        datadf = estimate_heterozygosity_from_af(datafile, windowsize, windowstep)
    elif stat == "PI":
        datadf = load_pi_data(datafile)
    else:
        print("oops")

    #saving plot
    plot_stat_w_groundtruth_inv(datadf, stat, chromosome, inv_map, windowsize)




#running program
main_func("allele_freq_chr1.frq", "Heterozygosity", "chr1", 1000, 100)
main_func("allele_freq_chr4.frq", "Heterozygosity", "chr4", 100000, 10000)
main_func("allele_freq_chr6.frq", "Heterozygosity", "chr6", 1000, 100)
main_func("allele_freq_chr14.frq", "Heterozygosity", "chr14", 1000, 100)
main_func("allele_freq_chr14.frq", "Heterozygosity", "chr14", 10000, 1000)
main_func("allele_freq_chr15.frq", "Heterozygosity", "chr15", 100000, 10000)


main_func("pi_win1000_chr1.windowed.pi", "PI", "chr1", 1000, 100)
main_func("pi_win100000_chr4.windowed.pi", "PI", "chr4", 100000, 10000)
main_func("pi_win1000_chr6.windowed.pi", "PI", "chr6", 1000, 100)
main_func("pi_win1000_chr14.windowed.pi", "PI", "chr14", 1000, 100)
main_func("pi_win10000_chr14.windowed.pi", "PI", "chr14", 10000, 1000)
main_func("pi_win100000_chr15.windowed.pi", "PI", "chr15", 100000, 10000)

