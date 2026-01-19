import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler



def _get_sitehaplist(record):
    """
    Helper function that take in a line (record) from vcf and Returns a site haplotype list needed for appending to matrix
    """

    GT_ls = record.strip().split("\t")[9:]
    hap_ls = []
    for gt in GT_ls:
        haps = gt.split("|")
        hap_ls.append(int(haps[0]))
        hap_ls.append(int(haps[1]))

    return hap_ls


def _get_record_position(record):
    """
    Helper function grabs the SNP coordinate for the input record and returns it as an integer
    """

    record_pos = int(record.strip().split("\t")[1])

    return record_pos



def pca_on_window(haplo_matrix):
    """
    Takes in a ndarray of haplotypes for a window and performs a principal component analysis.
    Input look like:
            site1, site2, ... siteN
    hap1
    hap2
    ...
    hapN

    Returns PCA 2D matrix and percent variance for each PC
    """

    #normalize
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(haplo_matrix)

    #PCA
    pca = PCA(n_components=2)
    PCs = pca.fit_transform(scaled_data)
    pcv = pca.explained_variance_ratio_

    #convert to df
    df = pd.DataFrame(PCs, columns=["PC1", "PC2"])


    return df, pcv[0], pcv[1]



def run_pca_on_SNP_windows(vcfgz, SNPwindow, output_dir=""):
    """
    Runs principal component analysis on haplotypes from a phased vcf in non-overlapping SNP windows

    The function will slide a window of SNPwindow size through the vcf, extract haplotypes, and run
    PCA of 2 PCs using sklearn. It will then output a csv file of PCs for each window will be saved
    to the given output_dir (default is the current dir). After the windowing complete, the function
    will also save a csv of explained variance of each PC for each window.
    """

    #empty lists for variance explained df
    win = []
    pos = []
    pc1v = []
    pc2v = []

    #open file
    with gzip.open(vcfgz, "rt") as vcf:

        #setting counters
        window_counter = 0

        SNP_counter = 0
        start_pos = 0
        end_pos = 0

        #empty matrix
        matrix = []

        #loop through lines
        for line in vcf:

            #skip header
            if line.startswith("#"):
                continue

            #run through SNPs
            else:
                
                #update counters
                if SNP_counter == 0:
                    start_pos = _get_record_position(line)
                SNP_counter += 1
                
                #adding record to hap matrix
                site_hap_ls = _get_sitehaplist(line)
                matrix.append(site_hap_ls)

                #if end of window is reached
                if SNP_counter == SNPwindow:
                    #logging meta data
                    win.append(window_counter)
                    end_pos = _get_record_position(line)
                    pos.append((start_pos + end_pos) // 2)

                    #converting matrix to array and transposing
                    matrix = np.array(matrix)
                    matrix = np.transpose(matrix)

                    #running PCA
                    pcdf, v1, v2 = pca_on_window(matrix)
                    pc1v.append(v1)
                    pc2v.append(v2)
                    pcdf.to_csv(f"{output_dir}win{window_counter}_windowsize{SNPwindow}_from_{start_pos}_to_{end_pos}_{vcfgz}".replace(".vcf.gz", ".csv"), index=False)

                    #resetting counts
                    SNP_counter = 0
                    start_pos = 0
                    end_pos = 0

                    #counting window
                    window_counter += 1

                    #clear matrix
                    matrix = []

                else:
                    continue

        #cleanup
        if SNP_counter != 0:
            win.append(window_counter)
            end_pos = _get_record_position(line)
            pos.append((start_pos + end_pos) // 2)

            matrix = np.array(matrix)
            matrix = np.transpose(matrix)

            pcdf, v1, v2 = pca_on_window(matrix)
            pc1v.append(v1)
            pc2v.append(v2)
            pcdf.to_csv(f"{output_dir}win{window_counter}_windowsize{SNPwindow}_from_{start_pos}_to_{end_pos}_{vcfgz}".replace(".vcf.gz", ".csv"), index=False)

            #clear matrix
            matrix = []

    #saving variance explained csv
    vdf = pd.DataFrame({"Window":win, "POS":pos, "PC1_variance":pc1v, "PC2_variance":pc2v})
    vdf.to_csv(f"{output_dir}var_explained_windowsize{SNPwindow}_{vcfgz}".replace(".vcf.gz", ".csv"), index=False)

    print(f"Done sliding window PCA for {vcfgz} with window size {SNPwindow}")










###25 SNP windows
#inversions
run_pca_on_SNP_windows("maf0.01_chr9_1MB_inv_1MB_2014048bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr9/inv/")
run_pca_on_SNP_windows("maf0.01_chr19_1MB_inv_1MB_2415129bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr19/inv/")
run_pca_on_SNP_windows("maf0.01_chr17_1MB_inv_1MB_3600721bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr17/inv/")
run_pca_on_SNP_windows("maf0.01_chr15_1MB_inv_1MB_4529599bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr15.1/inv/")
run_pca_on_SNP_windows("maf0.01_chr15_1MB_inv_1MB_8082027bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr15.2/inv/")

#controls
run_pca_on_SNP_windows("maf0.01_chr9_control1_2014048bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr9/control1/")
run_pca_on_SNP_windows("maf0.01_chr9_control2_2014048bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr9/control2/")
run_pca_on_SNP_windows("maf0.01_chr19_control1_2415129bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr19/control1/")
run_pca_on_SNP_windows("maf0.01_chr19_control2_2415129bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr19/control2/")
run_pca_on_SNP_windows("maf0.01_chr17_control1_3600721bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr17/control1/")
run_pca_on_SNP_windows("maf0.01_chr17_control2_3600721bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr17/control2/")
run_pca_on_SNP_windows("maf0.01_chr15_control1_4529599bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr15.1/control1/")
run_pca_on_SNP_windows("maf0.01_chr15_control2_4529599bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr15.1/control2/")
run_pca_on_SNP_windows("maf0.01_chr15_control1_8082027bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr15.2/control1/")
run_pca_on_SNP_windows("maf0.01_chr15_control2_8082027bp.vcf.gz", 25, output_dir="pca_output/maf0.01_filtered/25SNP_windows/chr15.2/control2/")




###50 SNP windows
#inversions
run_pca_on_SNP_windows("maf0.01_chr9_1MB_inv_1MB_2014048bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr9/inv/")
run_pca_on_SNP_windows("maf0.01_chr19_1MB_inv_1MB_2415129bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr19/inv/")
run_pca_on_SNP_windows("maf0.01_chr17_1MB_inv_1MB_3600721bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr17/inv/")
run_pca_on_SNP_windows("maf0.01_chr15_1MB_inv_1MB_4529599bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr15.1/inv/")
run_pca_on_SNP_windows("maf0.01_chr15_1MB_inv_1MB_8082027bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr15.2/inv/")

#controls
run_pca_on_SNP_windows("maf0.01_chr9_control1_2014048bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr9/control1/")
run_pca_on_SNP_windows("maf0.01_chr9_control2_2014048bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr9/control2/")
run_pca_on_SNP_windows("maf0.01_chr19_control1_2415129bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr19/control1/")
run_pca_on_SNP_windows("maf0.01_chr19_control2_2415129bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr19/control2/")
run_pca_on_SNP_windows("maf0.01_chr17_control1_3600721bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr17/control1/")
run_pca_on_SNP_windows("maf0.01_chr17_control2_3600721bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr17/control2/")
run_pca_on_SNP_windows("maf0.01_chr15_control1_4529599bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr15.1/control1/")
run_pca_on_SNP_windows("maf0.01_chr15_control2_4529599bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr15.1/control2/")
run_pca_on_SNP_windows("maf0.01_chr15_control1_8082027bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr15.2/control1/")
run_pca_on_SNP_windows("maf0.01_chr15_control2_8082027bp.vcf.gz", 50, output_dir="pca_output/maf0.01_filtered/50SNP_windows/chr15.2/control2/")


###500 SNP windows
#inversions
run_pca_on_SNP_windows("maf0.01_chr9_1MB_inv_1MB_2014048bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr9/inv/")
run_pca_on_SNP_windows("maf0.01_chr19_1MB_inv_1MB_2415129bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr19/inv/")
run_pca_on_SNP_windows("maf0.01_chr17_1MB_inv_1MB_3600721bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr17/inv/")
run_pca_on_SNP_windows("maf0.01_chr15_1MB_inv_1MB_4529599bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr15.1/inv/")
run_pca_on_SNP_windows("maf0.01_chr15_1MB_inv_1MB_8082027bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr15.2/inv/")

#controls
run_pca_on_SNP_windows("maf0.01_chr9_control1_2014048bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr9/control1/")
run_pca_on_SNP_windows("maf0.01_chr9_control2_2014048bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr9/control2/")
run_pca_on_SNP_windows("maf0.01_chr19_control1_2415129bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr19/control1/")
run_pca_on_SNP_windows("maf0.01_chr19_control2_2415129bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr19/control2/")
run_pca_on_SNP_windows("maf0.01_chr17_control1_3600721bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr17/control1/")
run_pca_on_SNP_windows("maf0.01_chr17_control2_3600721bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr17/control2/")
run_pca_on_SNP_windows("maf0.01_chr15_control1_4529599bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr15.1/control1/")
run_pca_on_SNP_windows("maf0.01_chr15_control2_4529599bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr15.1/control2/")
run_pca_on_SNP_windows("maf0.01_chr15_control1_8082027bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr15.2/control1/")
run_pca_on_SNP_windows("maf0.01_chr15_control2_8082027bp.vcf.gz", 500, output_dir="pca_output/maf0.01_filtered/500SNP_windows/chr15.2/control2/")





























""" 
###200 SNP windows
#inversions
run_pca_on_SNP_windows("chr9_1MB_inv_1MB_2014048bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr9/inv/")
run_pca_on_SNP_windows("chr19_1MB_inv_1MB_2415129bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr19/inv/")
run_pca_on_SNP_windows("chr17_1MB_inv_1MB_3600721bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr17/inv/")
run_pca_on_SNP_windows("chr15_1MB_inv_1MB_4529599bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr15/inv1/")
run_pca_on_SNP_windows("chr15_1MB_inv_1MB_8082027bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr15/inv2/")

#controls
run_pca_on_SNP_windows("chr9_control1_2014048bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr9/control1/")
run_pca_on_SNP_windows("chr9_control2_2014048bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr9/control2/")
run_pca_on_SNP_windows("chr19_control1_2415129bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr19/control1/")
run_pca_on_SNP_windows("chr19_control2_2415129bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr19/control2/")
run_pca_on_SNP_windows("chr17_control1_3600721bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr17/control1/")
run_pca_on_SNP_windows("chr17_control2_3600721bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr17/control2/")
run_pca_on_SNP_windows("chr15_control1_4529599bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr15/control1.1/")
run_pca_on_SNP_windows("chr15_control2_4529599bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr15/control2.1/")
run_pca_on_SNP_windows("chr15_control1_8082027bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr15/control1.2/")
run_pca_on_SNP_windows("chr15_control2_8082027bp.vcf.gz", 200, output_dir="pca_output/200SNP_windows/chr15/control2.2/")




###400 SNP windows
#inversions
run_pca_on_SNP_windows("chr9_1MB_inv_1MB_2014048bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr9/inv/")
run_pca_on_SNP_windows("chr19_1MB_inv_1MB_2415129bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr19/inv/")
run_pca_on_SNP_windows("chr17_1MB_inv_1MB_3600721bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr17/inv/")
run_pca_on_SNP_windows("chr15_1MB_inv_1MB_4529599bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr15/inv1/")
run_pca_on_SNP_windows("chr15_1MB_inv_1MB_8082027bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr15/inv2/")

#controls
run_pca_on_SNP_windows("chr9_control1_2014048bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr9/control1/")
run_pca_on_SNP_windows("chr9_control2_2014048bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr9/control2/")
run_pca_on_SNP_windows("chr19_control1_2415129bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr19/control1/")
run_pca_on_SNP_windows("chr19_control2_2415129bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr19/control2/")
run_pca_on_SNP_windows("chr17_control1_3600721bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr17/control1/")
run_pca_on_SNP_windows("chr17_control2_3600721bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr17/control2/")
run_pca_on_SNP_windows("chr15_control1_4529599bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr15/control1.1/")
run_pca_on_SNP_windows("chr15_control2_4529599bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr15/control2.1/")
run_pca_on_SNP_windows("chr15_control1_8082027bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr15/control1.2/")
run_pca_on_SNP_windows("chr15_control2_8082027bp.vcf.gz", 400, output_dir="pca_output/400SNP_windows/chr15/control2.2/")


###4000 SNP windows
#inversions
run_pca_on_SNP_windows("chr9_1MB_inv_1MB_2014048bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr9/inv/")
run_pca_on_SNP_windows("chr19_1MB_inv_1MB_2415129bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr19/inv/")
run_pca_on_SNP_windows("chr17_1MB_inv_1MB_3600721bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr17/inv/")
run_pca_on_SNP_windows("chr15_1MB_inv_1MB_4529599bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr15/inv1/")
run_pca_on_SNP_windows("chr15_1MB_inv_1MB_8082027bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr15/inv2/")

#controls
run_pca_on_SNP_windows("chr9_control1_2014048bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr9/control1/")
run_pca_on_SNP_windows("chr9_control2_2014048bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr9/control2/")
run_pca_on_SNP_windows("chr19_control1_2415129bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr19/control1/")
run_pca_on_SNP_windows("chr19_control2_2415129bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr19/control2/")
run_pca_on_SNP_windows("chr17_control1_3600721bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr17/control1/")
run_pca_on_SNP_windows("chr17_control2_3600721bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr17/control2/")
run_pca_on_SNP_windows("chr15_control1_4529599bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr15/control1.1/")
run_pca_on_SNP_windows("chr15_control2_4529599bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr15/control2.1/")
run_pca_on_SNP_windows("chr15_control1_8082027bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr15/control1.2/")
run_pca_on_SNP_windows("chr15_control2_8082027bp.vcf.gz", 4000, output_dir="pca_output/4000SNP_windows/chr15/control2.2/")
"""