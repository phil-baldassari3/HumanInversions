import pandas as pd
import numpy as np



def remove_bad_regions(bedgraph_df, region_list):
    """
    Function takes in the DataFrame generated from compute_repeat_and_SNP_density and removes specified regions.
    This is used to remove centromeric regions that skew the branch lengths.
    
    :param bedgraph_df: DataFrame from compute_repeat_and_SNP_density
    :param region_list: List of tuples each containing 3 elements: "chrom", start_site, end_site

    :type df: DataFrame
    :type region_list: List of Tuples

    :return: DataFrame with removed regions
    :rtype: DataFrame
    """

    #setting empty mask
    keep_rows = np.ones(len(bedgraph_df), dtype=bool)

    #masking rows for windows entirely or partially overlapping with regions to remove
    for chrom, region_start, region_end in region_list:
        keep_rows &= ~((bedgraph_df["CHROM"] == chrom) & (bedgraph_df["START"] < region_end) & (bedgraph_df["END"] > region_start))

    #removing bad rows
    new_df = bedgraph_df[keep_rows]


    return new_df




def main(bedgraph, remove_regions, outfile):
    """
    Main function that runs remove_bad_regions and outputs a bedgraph file
    
    :param bedfile: bedfile name
    :type bedfile: str
    :param remove_regions: list of tuples with regions to remove
    :type remove_regions: List of Tuples
    :param outfile: name of output TSV
    :type outfile: str

    :return: None
    """

    #opening file
    df = pd.read_csv(bedgraph, sep="\t", header=None, names=["CHROM", "START", "END", "PARAMETER"])

    cleaned_df = remove_bad_regions(df, remove_regions)

    cleaned_df.to_csv(outfile, sep="\t", index=False, header=False)



### REGIONS TO REMOVE ###
#centromere coordinates from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/
#removing pericentromere 3.5Mb on either side of centromere or 1Mb on either side of the centromere
#also removing large gap on Y
remove_these35 = [
    ("chr1", 122026460-3500000, 125184587+3500000),
    ("chr2", 92188146-3500000, 94090557+3500000),
    ("chr3", 90772459-3500000, 93655574+3500000),
    ("chr4", 49708101-3500000, 51743951+3500000),
    ("chr5", 46485901-3500000, 50059807+3500000),
    ("chr6", 58553889-3500000, 59829934+3500000),
    ("chr7", 58169654-3500000, 60828234+3500000),
    ("chr8", 44033745-3500000, 45877265+3500000),
    ("chr9", 43236168-3500000, 45518558+3500000),
    ("chr10", 39686683-3500000, 41593521+3500000),
    ("chr11", 51078349-3500000, 54425074+3500000),
    ("chr12", 34769408-3500000, 37185252+3500000),
    ("chr13", 16000001-3500000, 18051248+3500000),
    ("chr14", 16000001-3500000, 18173523+3500000),
    ("chr15", 17000001-3500000, 19725254+3500000),
    ("chr16", 36311159-3500000, 38280682+3500000),
    ("chr17", 22813680-3500000, 26885980+3500000),
    ("chr18", 15460900-3500000, 20861206+3500000),
    ("chr19", 24498981-3500000, 27190874+3500000),
    ("chr20", 26436233-3500000, 30038348+3500000),
    ("chr21", 10864561-3500000, 12915808+3500000),
    ("chr22", 12954789-3500000, 15054318+3500000),
    ("chrX", 58605580-3500000, 62412542+3500000),
    ("chrY", 10316945-3500000, 10544039+3500000),
    ("chrY", 26784117, 56558617)
]


remove_these1 = [
    ("chr1", 122026460-1000000, 125184587+1000000),
    ("chr2", 92188146-1000000, 94090557+1000000),
    ("chr3", 90772459-1000000, 93655574+1000000),
    ("chr4", 49708101-1000000, 51743951+1000000),
    ("chr5", 46485901-1000000, 50059807+1000000),
    ("chr6", 58553889-1000000, 59829934+1000000),
    ("chr7", 58169654-1000000, 60828234+1000000),
    ("chr8", 44033745-1000000, 45877265+1000000),
    ("chr9", 43236168-1000000, 45518558+1000000),
    ("chr10", 39686683-1000000, 41593521+1000000),
    ("chr11", 51078349-1000000, 54425074+1000000),
    ("chr12", 34769408-1000000, 37185252+1000000),
    ("chr13", 16000001-1000000, 18051248+1000000),
    ("chr14", 16000001-1000000, 18173523+1000000),
    ("chr15", 17000001-1000000, 19725254+1000000),
    ("chr16", 36311159-1000000, 38280682+1000000),
    ("chr17", 22813680-1000000, 26885980+1000000),
    ("chr18", 15460900-1000000, 20861206+1000000),
    ("chr19", 24498981-1000000, 27190874+1000000),
    ("chr20", 26436233-1000000, 30038348+1000000),
    ("chr21", 10864561-1000000, 12915808+1000000),
    ("chr22", 12954789-1000000, 15054318+1000000),
    ("chrX", 58605580-1000000, 62412542+1000000),
    ("chrY", 10316945-1000000, 10544039+1000000),
    ("chrY", 26784117, 56558617)
]





main("All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph", remove_these1, "Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount.bedgraph")
main("All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph", remove_these1, "Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount.bedgraph")
main("ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph", remove_these1, "Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.bedgraph")
main("ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph", remove_these1, "Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.bedgraph")

