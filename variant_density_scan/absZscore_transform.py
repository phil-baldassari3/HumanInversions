import pandas as pd
import numpy as np


def _chrom_sorter(chrom_name):
    """
    Helper function used to sort chromosome names by number
    """

    chrom_num = int(chrom_name.replace("chr", ""))

    return chrom_num



def absZtransform_data(bedgraph, outbedgraph):
    """
    Function takes in a bedgraph file, Z-transforms the data column, and outputs a new bedgraph.
    Note that the output data are the absolute value Z-scores. The data is Z-transformed separately for each chromosome.

    :param bedgraph: input bedgraph file to be Z-transformed
    :param outbedgraph: output bedgraph base file name

    :type bedgraph: str
    :type outbedgraph: str
    """

    print(f"\nAbsolute Value Z-Score Transform on {bedgraph}")

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
    chrom_df_list = []
    for chrom in chroms:
        #filtering df
        filtered_df = df[df["CHROM"] == chrom].copy()
        
        #computing mean and sd
        m = filtered_df["Data"].mean()
        s = filtered_df["Data"].std()
        print(f"{chrom}: mean={m}, sd={s}")

        #Z-transform
        filtered_df["Z"] = np.abs((filtered_df["Data"] - m) / s)

        #append df to list
        chrom_df_list.append(filtered_df)

    #concatenate dfs
    final_df = pd.concat(chrom_df_list)
    final_df = final_df[["CHROM", "START", "END", "Z"]]

    #outputting file
    final_df.to_csv(f"{outbedgraph}.bedgraph", index=False, header=False, sep="\t")




absZtransform_data("win500Kb_step50Kb_0.01_all_biallelic_variant_count.bedgraph", "win500Kb_step50Kb_0.01_all_biallelic_variant_Z")
absZtransform_data("win100Kb_step10Kb_0.01_all_biallelic_variant_count.bedgraph", "win100Kb_step10Kb_0.01_all_biallelic_variant_Z")
absZtransform_data("win10Kb_step1Kb_0.01_all_biallelic_variant_count.bedgraph", "win10Kb_step1Kb_0.01_all_biallelic_variant_Z")


absZtransform_data("win500Kb_step50Kb_0.01_SNP_count.bedgraph", "win500Kb_step50Kb_0.01_SNP_Z")
absZtransform_data("win100Kb_step10Kb_0.01_SNP_count.bedgraph", "win100Kb_step10Kb_0.01_SNP_Z")
absZtransform_data("win10Kb_step1Kb_0.01_SNP_count.bedgraph", "win10Kb_step1Kb_0.01_SNP_Z")


absZtransform_data("win500Kb_step50Kb_0.01_SNV_count.bedgraph", "win500Kb_step50Kb_0.01_SNV_Z")
absZtransform_data("win100Kb_step10Kb_0.01_SNV_count.bedgraph", "win100Kb_step10Kb_0.01_SNV_Z")
absZtransform_data("win10Kb_step1Kb_0.01_SNV_count.bedgraph", "win10Kb_step1Kb_0.01_SNV_Z")