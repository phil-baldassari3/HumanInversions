import pandas as pd
import numpy as np


def _window_size_bookkeeping(win, step, mult):
    """
    Helper function to keep track of how many input windows are in each varwindow and varstep.

    :param win: size of windows in bp in the input bedgraph
    :param step: size of steps in bp in the input bedgraph
    :param mult: how many times larger should the variance windows be

    :type win: int
    :type step: int
    :type mult: int

    :returns: number of windows in varwindow, number of windows in varstep
    :rtype: int, int
    """

    #finding var win and step sizes in bp
    varwin = win * mult
    varstep = step * mult

    #finding number of windows in varwin
    numwin_in_varwin = varwin // step

    #finding number of windows in varstep
    numwin_in_varstep = varstep // step

    return numwin_in_varwin, numwin_in_varstep





def _chrom_sorter(chrom_name):
    """
    Helper function used to sort chromosome names by number
    """

    chrom_num = int(chrom_name.replace("chr", ""))

    return chrom_num



def var_scan(bedgraph, input_window_size, input_step_size, window_multiplier, outbedgraph):
    """
    Function takes in a bedgraph file of windowed data sorted by chromosome and position and with equal sized windows
    and computes the variance of the data in a window size and step a defined magnitude greater than the original window size

    :param bedgraph: input bedgraph file to variance scanned
    :param input_window_size: size of windows in bp in the input bedgraph
    :param input_step_size: size of steps in bp in the input bedgraph
    :param window_multiplier: how many times larger should the variance windows be
    :param outbedgraph: output bedgraph base file name

    :type bedgraph: str
    :type input_window_size: int
    :type input_step_size: int
    :type window_multiplier: int
    :type outbedgraph: str
    """

    print(f"\nVariance scan on {bedgraph} in {input_window_size*window_multiplier} bp windows")

    #finding number of windows in each varwin and varstep
    varwin_wins, varstep_wins = _window_size_bookkeeping(input_window_size, input_step_size, window_multiplier)
    print(f"Expecting {varwin_wins} rows per variance window ({varwin_wins * input_step_size:,} bp)")

    #opening bedgraph as a dataframe
    df = pd.read_csv(bedgraph, sep="\t", header=None, names=["CHROM", "START", "END", "Data"])
    df.dropna(inplace=True)

    #find chromosomes
    temp_chroms = set(df["CHROM"])

    #sort chromosomes
    autosomes = [c for c in temp_chroms if c not in ("chrX", "chrY")]
    autosomes.sort(key=_chrom_sorter)
    chroms = autosomes
    if "chrX" in temp_chroms:
        chroms.append("chrX")
    if "chrY" in temp_chroms:
        chroms.append("chrY")

    #iterating through chromosomes and outputting chromosome dfs with Z-transformed data
    newchrom = []
    newstart = []
    newend = []
    vars = []
    for chrom in chroms:
        #filtering df
        filtered_df = df[df["CHROM"] == chrom].reset_index(drop=True)
        print(chrom)
        #windowing
        for startidx in range(0, len(filtered_df)-varwin_wins, varstep_wins):
            chunk = filtered_df.iloc[startidx:startidx+varwin_wins]

            first_row = chunk.iloc[0]
            last_row  = chunk.iloc[-1]

            newchrom.append(first_row["CHROM"])
            newstart.append(first_row["START"])
            newend.append(last_row["END"])
            vars.append(np.var(chunk["Data"], ddof=1))


    #outputting file
    var_df = pd.DataFrame({"CHROM":newchrom, "START":newstart, "END":newend, "Var":vars})
    var_df.to_csv(f"{outbedgraph}.bedgraph", index=False, header=False, sep="\t")




var_scan("win100Kb_step10Kb_0.01_all_biallelic_variant_count.bedgraph", 100000, 10000, 10, "varwin1Mb_win100Kb_step10Kb_0.01_all_biallelic_variant_Var")
var_scan("win10Kb_step1Kb_0.01_all_biallelic_variant_count.bedgraph", 10000, 1000, 10, "varwin100Kb_win10Kb_step1Kb_0.01_all_biallelic_variant_Var")


var_scan("win100Kb_step10Kb_0.01_SNP_count.bedgraph", 100000, 10000, 10, "varwin1Mb_win100Kb_step10Kb_0.01_SNP_Var")
var_scan("win10Kb_step1Kb_0.01_SNP_count.bedgraph", 10000, 1000, 10, "varwin100Kb_win10Kb_step1Kb_0.01_SNP_Var")


var_scan("win100Kb_step10Kb_0.01_SNV_count.bedgraph", 100000, 10000, 10, "varwin1Mb_win100Kb_step10Kb_0.01_SNV_Var")
var_scan("win10Kb_step1Kb_0.01_SNV_count.bedgraph", 10000, 1000, 10, "varwin100Kb_win10Kb_step1Kb_0.01_SNV_Var")