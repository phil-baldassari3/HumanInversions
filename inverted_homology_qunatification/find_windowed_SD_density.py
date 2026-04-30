import os
import pandas as pd
import subprocess


def _make_bedfile_of_windows(chrom_lengths, window_size, window_step):
    """
    Helper function saves a bedfile to the working directory that windows the genome with a 
    specified window size and step. To do this, the function needs a text file where each line
    contains a chromosome name and its length separated by a tab: e.g. chr1\t248956422

    :param : file name for tab separated table of chromosome sizes
    :param : window size in bp
    :param : window step in bp

    :type : str
    :type : int
    :type : int
    """

    #open files
    with open("window.bed", "w") as winbed, open(chrom_lengths, "r") as chromlens:

        #iterate through chromosomes and grab chrom and length
        for line in chromlens:
            linelist = line.strip().split("\t")
            chrom = linelist[0]
            length = int(linelist[1])

            #iterate through windows and write to bedfile
            for start in range(0, length, window_step):
                end = start + window_size
                bedstr = f"{chrom}\t{start}\t{end}\n"

                winbed.write(bedstr)


def _count_intersecting_SD_bps(windowbed, mergedSDbed):
    """
    Helper function takes a bedfile of windows and a bedfile of merged Segmental Duplications
    and counts the number of SD bps in each window using bedtools intersect -wao run as a 
    subprocess. An output bedfile of intersects and bp counts is saved to the current directory.

    :param windowbed: bedfile name of windows across the genome made from `_make_bedfile_of_windows`
    :param mergedSDbed: bedfile name of merged SD regions

    :type windowbed: str
    :type mergedSDbed: str
    """

    #write intersect bedfile
    with open("windowed_SD_intersects.bed", "w") as outbed:
        subprocess.Popen(["bedtools", "intersect", "-a", windowbed, "-b", mergedSDbed, "-wao"], stdout=outbed).wait()


def _compute_windowed_SD_density(intersectbed, window_size):
    """
    Helper function takes the intersect bedfile output from `_count_intersecting_SD_bps`,
    converts the bed to a pandas dataframe and computes the SD density per window as a 
    proportion of SD_bp_count/window_size.

    :param intersectbed: bedfile name of file output of `_count_intersecting_SD_bps`
    :param window_size: window size in bp

    :type intersectbed: str
    :type window_size: int

    :returns: Dataframe of windowed SD density
    :rtype: pandas.DataFrame
    """

    #open file as df and select necessary columns
    df = pd.read_csv(intersectbed, sep="\t", names=["CHROM", "START", "END", "CHROM2", "START2", "END2", "SD_BP"])
    df = df[["CHROM", "START", "END", "SD_BP"]]

    #group spans that intersect with multiple SD regions and sum the BP count
    df = df.groupby(by=["CHROM", "START", "END"], as_index=False).sum()

    #compute density
    df["SD_density"] = df["SD_BP"] / window_size

    #select necessary columns
    df = df[["CHROM", "START", "END", "SD_density"]]

    return df



def run_the_program(chrom_sizes, mergedSDs, winsize, winstep, output):
    """
    Main Function that takes in a file of chromosome sizes, makes a temporary bedfile of
    sliding windows defined by a size and step, counts the number of intersecting SD bps
    in each window, and computes the proportion of each window made up of SDs. The function
    outputs a bedgraph file with results and deletes the temporary files it made along the way.

    :param chrom_sizes: file name for tab separated table of chromosome sizes, see `_make_bedfile_of_windows`
    :param mergedSDs: bedfile name of merged SD regions
    :param winsize: window size in bp
    :param winstep: window step in bp
    :param output: base name for output bedgraph

    :type chrom_sizes: str
    :type mergedSDs: str
    :type winsize: int
    :type winstep: int
    :type output: str
    """

    #make windows, makes "window.bed" file
    _make_bedfile_of_windows(chrom_sizes, winsize, winstep)

    #count intersecting bps, makes "windowed_SD_intersects.bed" file
    _count_intersecting_SD_bps("window.bed", mergedSDs)

    #compute density 
    den_df = _compute_windowed_SD_density("windowed_SD_intersects.bed", winsize)

    #delete temporary files
    if os.path.exists("window.bed"):
        os.remove("window.bed")
    if os.path.exists("windowed_SD_intersects.bed"):
        os.remove("windowed_SD_intersects.bed")

    #output bedgraph
    den_df.to_csv(f"{output}.bedgraph", sep="\t", header=False, index=False)





def main():
    run_the_program("hg38_chrom_1_to_22_X_and_Y_sizes.txt", "vollger_hg38_SDs.merged.bed", 100000, 10000, "windowed_SD_density")

if __name__ == "__main__":
    main()