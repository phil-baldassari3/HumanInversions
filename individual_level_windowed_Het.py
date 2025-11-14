import gzip
import csv
import numpy as np
from statistics import geometric_mean
from scipy.stats import pmean
from multiprocessing import Pool


def vcf2persiteHet(vcfgz):
    """
    Loads in a gzipped vcf and returns a list of lists.
    The data is 0 for a homozygous site and 1 for a 
    heterozygous site.

    Returns
    chromosome str
    header list
    ###output file name
    np array: [[pos, het1, het2...],]
    """

    print(f"Loading data from {vcfgz}...")

    #Hom Het matrix
    homhet_data = []

    #output persite heterozygosity file
    #with open(persitefile, "a") as hetfile:
    #opening file and reading line by line
    with gzip.open(vcfgz, "rt") as vcf:
        for line in vcf:

            if line.startswith("##"):
                continue
            elif line.startswith("#CHROM"):
                editline = line.replace("#", "").rstrip()
                line2list = editline.split("\t")
                
                header = line2list[0:2] + line2list[9:]
                header += ["Geometric_mean", "p10_mean"]
            
            else:
                rowlist = []
                line2list = line.split("\t")
                #strip end of the line
                line2list[-1] = line2list[-1].strip()
                #grab chromosome
                chrom = (line2list[0])
                #constructing row
                rowlist.append(int(line2list[1]))
                for gt in line2list[9:]:
                    if gt == "0|0" or gt == "1|1":
                        rowlist.append(0)
                    elif gt == "0|1" or gt == "1|0":
                        rowlist.append(1)

                #write row to file
                #writer = csv.writer(hetfile)
                #writer.writerow(rowlist)
                    
                homhet_data.append(rowlist)

    #convert to array
    homhet_data = np.array(homhet_data)

    print(f"FINISHED loading data from {vcfgz}")

    return chrom, header, homhet_data   #persitefile




def persite2windowedHet(matrix, chromosome, header, SNPwindow, SNPstep):    #persitefile
    """
    Function takes the HomHet data matrix and returns windowed heterozygosity
    for each sample and the geometric mean of heterozygosity across samples.
    Windows are SNP-based

    appends to output file: chr, pos, indv_windowed_het...geomean, p10mean
    """

    print(f"{SNPwindow} SNP sliding window with {SNPstep} for {chromosome}...")

    #creating output file
    outfilename = f"Heterozygosity_{chromosome}_{SNPwindow}_SNP_windows.csv"
    with open(outfilename, "w") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)

    #opening array from file
    #matrix = np.loadtxt(persitefile, delimiter=',')
    #sliding window
    for start in range(0, len(matrix), SNPstep):

        #setting windowing
        end = start + SNPwindow

        #window chunk
        chunk = matrix[start:end, :]
        startpos = chunk[0, 0]
        endpos = chunk[-1, 0]
        chunkpos = (startpos + endpos) // 2
        chunk_length = (endpos - startpos) + 1

        #summing per site heterozygosity
        winHets = np.sum(chunk[:, 1:], axis=0, dtype=np.float64)
        
        #normalizing by window length
        winHets /= chunk_length
        #compute geometrix mean
        geomean = geometric_mean([x for x in winHets if x > 0])
        #compute p^10 mean
        p10mean = pmean(winHets, 10)

        #outputting to file
        output_row = [chromosome, chunkpos] + list(winHets) + [geomean, p10mean]

        with open(outfilename, "a") as outf:
            writer = csv.writer(outf)
            writer.writerow(output_row)

    print(f"FINISHED {SNPwindow} SNP sliding window for {chromosome}.")





def main_func(vcffile):
    """Runs the program with 3 diffrerent window sizes: 100, 1000, 10,000"""

    chrom, headerline, data = vcf2persiteHet(vcffile)
    persite2windowedHet(data, chrom, headerline, 100, 10)
    persite2windowedHet(data, chrom, headerline, 1000, 100)
    persite2windowedHet(data, chrom, headerline, 10000, 1000)


#main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr4.vcf.gz")
#main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr6.vcf.gz")
#main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz")
#main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz")
#main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr10.vcf.gz")


main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz")
main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz")
main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz")
main_func("beagle_phased_biallelic_SNPs_1000GP30X_chr14.vcf.gz")


"""
#for parallel mapping
vcfs = [
    "beagle_phased_biallelic_SNPs_1000GP30X_chr4.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr6.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr8.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr10.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr14.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz",
    "beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz",
]

#run in parallel
def run_in_parallel():
    pool = Pool(processes=2)
    pool.map(main_func, vcfs)


if __name__ == '__main__':
    run_in_parallel()






chrom, headerline, persitefilename = vcf2persiteHet("5samples_100000SNPs_beagle_benchmark.vcf.gz")
#chrom, outputf, data = vcf2persiteHet("beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz")

persite2windowedHet(persitefilename, chrom, headerline, 1000, 100)

"""
