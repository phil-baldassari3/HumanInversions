import pandas as pd



def _load_bedgraph(bedgraph, dataname):
    """
    Helper function that loads in a bedgraph as a dataframe.
    """

    df = pd.read_csv(bedgraph, sep="\t", header=None, names=["CHROM", "START", "END", dataname])

    return df


def compute_SNP_to_SNV_ratio(SNPbedgraph, SNVbedgraph):

    #open bedgraphs
    snps = _load_bedgraph(SNPbedgraph, "SNP_count")
    snvs = _load_bedgraph(SNVbedgraph, "SNV_count")

    #merge dataframes
    merged_df = pd.merge(snps, snvs, on=["CHROM", "START", "END"], how='inner')
    assert(len(merged_df) == len(snps))
    assert(len(merged_df) == len(snvs))

    #compute ratio
    merged_df = merged_df[merged_df["SNV_count"] != 0]
    merged_df["SNP_to_SNV_ratio"] = merged_df["SNP_count"] / merged_df["SNV_count"]

    #output new file
    merged_df = merged_df[["CHROM", "START", "END", "SNP_to_SNV_ratio"]]
    merged_df.to_csv(SNPbedgraph.replace("SNP_count", "SNPtoSNV_ratio"), sep="\t", index=False, header=False)



compute_SNP_to_SNV_ratio("win500Kb_step50Kb_0.01_SNP_count.bedgraph", "win500Kb_step50Kb_0.01_SNV_count.bedgraph")
compute_SNP_to_SNV_ratio("win100Kb_step10Kb_0.01_SNP_count.bedgraph", "win100Kb_step10Kb_0.01_SNV_count.bedgraph")
compute_SNP_to_SNV_ratio("win10Kb_step1Kb_0.01_SNP_count.bedgraph", "win10Kb_step1Kb_0.01_SNV_count.bedgraph")