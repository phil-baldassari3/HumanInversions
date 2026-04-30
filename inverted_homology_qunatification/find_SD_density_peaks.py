import os
import pandas as pd
import subprocess


def _filter_windowed_SD_density(density_bedgraph, threshold):
    """
    Helper function takes in the windowed SD density bedgraph output from `find_windowed_SD_density.py`
    and returns a dataframe of windowed data above a specified threshold (recommended 0.25)

    :param density_bedgraph: name of bedgraph file of windowed SD density
    :param threshold: SD density threshold to filter for windows above, recommended 0.25

    :type density_bedgraph: str
    :type threshold: float

    :returns: filtered dataframe
    :rtype: pandas.DataFrame
    """

    #open bedgraph
    df = pd.read_csv(density_bedgraph, sep="\t", names=["CHROM", "START", "END", "SD_density"])

    #filter by threshold
    df = df[df["SD_density"] >= threshold]

    return df


def _merge_filtered_windows(filtered_df, output_bed):
    """
    Helper function converts the filtered dataframe from `_filter_windowed_SD_density` to a
    bedfile of windows, then runs bedtools merge to merge adjacent windows together into peak
    spans, and saves an output bedgraph of high SD density peaks.

    :param filtered_df: dataframe returned from `_filter_windowed_SD_density`
    :param output_bed: base name for output bed file

    :type filtered_df: pandas.DataFrame
    :type output_bed: str
    """

    #grab bed only columns
    bed_df = filtered_df[["CHROM", "START", "END"]]

    #saving temporary bedfile
    bed_df.to_csv(f"TEMP_{output_bed}.bed", sep="\t", header=False, index=False)

    #merging windows into spans
    with open(f"{output_bed}.bed", "w") as outbed:
        subprocess.Popen(["bedtools", "merge", "-i", f"TEMP_{output_bed}.bed"], stdout=outbed).wait()

    #delete temporary file
    if os.path.exists(f"TEMP_{output_bed}.bed"):
        os.remove(f"TEMP_{output_bed}.bed")



def run_the_program(densitybedgraph, thresh, output):
    """
    Main function that takes in a windowed SD density bedgraph from `find_windowed_SD_density.py`,
    filters for windows above a specified density threshold, and outputs a bedfiles of filtered 
    windows merged into peak spans.

    :param densitybedgraph: name of bedgraph file of windowed SD density
    :param thresh: SD density threshold to filter for windows above, recommended 0.25
    :param output: base name for output bed file

    :type densitybedgraph: str
    :type thresh: float
    :type output: str
    """

    filtered_data = _filter_windowed_SD_density(densitybedgraph, thresh)
    _merge_filtered_windows(filtered_data, output)



def main():
    run_the_program("windowed_SD_density.bedgraph", 0.25, "SD_density_peaks")

if __name__ == "__main__":
    main()