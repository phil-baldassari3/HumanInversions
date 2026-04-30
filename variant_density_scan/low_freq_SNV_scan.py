"""
This script takes in a vcf.gz of biallelic SNPs and counts how many low frequency SNVs are 
in each window with a defined <=frequency, window size in bp and window step in bp. 
The program outputs a csv with MAF for each SNP in the vcf and a bedfile with the windowed results.
"""

import gzip
from collections import deque


### CLASSES ###
class MafCsvRowLoader():
    def __init__(self, csvrow):
        """
        This class is for parsing each row of the csv output from `persite_MAF` for it to be added to a deque.
        
        :param csvrow: row from the csv file output by `persite_MAF`

        :type csvrow: str

        :var chrom: row's chromosome
        :var pos: row's variant position
        :var maf: variant's minor allele frequency

        :vartype chrom: str
        :vartype pos: int
        :vartype maf: float
        """

        row_list = csvrow.strip().split(",")
        self.chrom = row_list[0]
        self.pos = int(row_list[1])
        self.maf = float(row_list[2])




## HELPER FUNCTIONS ##

def _joinany(sep, mylist):
    """
    Helper function to join a list by a sep regardless of element types
    """

    strs = [str(x) for x in mylist]

    return sep.join(strs)



def _compute_record_MAF(record_GTs):
    """
    Helper function to compute teh minor allele frequency for a record. The function takes
    a string of genotypes from the record. THis helper function is used in `persite_MAF`.

    :param record_GTs: string of record genotypes that look like: "0|0", "1/0"...
    :type record_GTs: str

    :returns: Minor allele frequency estimate for that record
    :rtype: float
    """

    #counting alleles
    num0s = record_GTs.count("0")
    num1s = record_GTs.count("1")
    total = num0s + num1s

    #computing maf
    maf = min(num0s, num1s) / total

    return maf


def _pop_bp_step(d, step):
    """
    Helper function used to pop a number of elements from the left
    of a deque equal to the step size in bp. This is done in place.
    """

    #count how many elements to pop
    count_ele_to_pop = 0
    for idx in range(len(d)):
        if d[idx].pos <= step:
            count_ele_to_pop += 1
        else:
            break
    
    #pop from left
    for i in range(count_ele_to_pop):
        d.popleft()



def _count_rare_SNVs(window_data, maf_thresh):
    """
    Helper function used to count the number of SNVs below the specified MAF threshold in a window.

    :param window_data: a window of data from csv parsed with `MafCsvRowLoader`
    :param maf_thresh: minor allele frequnecy threshold below which to count variants

    :type window_data: deque
    :type maf_thresh: float

    :returns: count of rare variants for that window
    :rtype: int
    """

    #iterating through window and counting rare variants
    rare_snv_count = 0
    for snv in window_data:
        if snv.maf <= maf_thresh and snv.maf > 0:
            rare_snv_count += 1

    return rare_snv_count


def _count_common_SNPs(window_data, maf_thresh):
    """
    Helper function used to count the number of SNPs above the specified MAF threshold in a window.

    :param window_data: a window of data from csv parsed with `MafCsvRowLoader`
    :param maf_thresh: minor allele frequnecy threshold below which to count variants

    :type window_data: deque
    :type maf_thresh: float

    :returns: count of common SNP variants for that window
    :rtype: int
    """

    #iterating through window and counting rare variants
    common_snp_count = 0
    for snv in window_data:
        if snv.maf >= maf_thresh:
            common_snp_count += 1

    return common_snp_count






### MAIN FUNCTIONS ###

def persite_MAF(vcfgz, GTindex=9):
    """
    Function takes a gzipped vcf file and estimates the minor allele frequency at each record.
    The output is a csv file with the same basename as the vcf + _MAF.csv saved to the same directory.
    The vcf MUST be ordered by chromosome, i.e. each chromosome's variants are together
    MAF = min(#of0s, #of1s) / (#of0s + #of1s)

    :param vcfgz: filename of input vcf.gz
    :param GTindex: optional argument to denote the index position of the starting Genotype column. Make sure you check this in your vcf and use python 0-based indexing
    
    :type vcfgz: str
    :type GTindex: int, optional
    """

    #opening files
    with gzip.open(vcfgz, "rt") as vcf, open(vcfgz.replace(".vcf.gz", "_MAF.csv"), "a") as outcsv:

        #write header line
        outcsv.write(",".join(["CHROM", "POS", "MAF"]) + "\n")

        #string for progress messaging
        chrom_prog = "reading in vcf file"

        #loop through lines
        for line in vcf:

            #skip header lines
            if line.startswith("#"):
                continue

            #record lines
            else:

                #grab data
                chrom = line.strip().split("\t")[0]
                pos = line.strip().split("\t")[1]
                gt_str = "\t".join(line.strip().split("\t")[GTindex:])

                #compute MAF
                record_maf = _compute_record_MAF(gt_str)

                #progess mesaging
                if chrom != chrom_prog:
                    print(f"Finished {chrom_prog}, Computing MAF for {chrom}")
                    chrom_prog = chrom

                #append to csv
                outcsv.write(_joinany(",", [chrom, pos, record_maf]) + "\n")




def sliding_window_SNV_count(MAFcsv, window_size, window_step, frequency, output):
    """
    Function performs a sliding window counting of variants below a specified MAF threshold
    using the csv file generated from `persite_MAF`. The function outputs a bedgraph file with results.
    Note that this output file has a header line starting with "#"

    :param MAFcsv: input minor allele frequency csv file
    :param window_size: size of sliding window in bp
    :param window_step: distance to slide window each iteration in bp
    :param frequency: minor allele frequency threshold below which to count the number of variants in the window
    :param output: basename for the output bedgraph file


    :type MAFcsv: str
    :type window_size: int
    :type window_step: int
    :type frequency: float
    :type output: str
    """

    print("Sliding window rare SNV count...")

    #nested functions that make the flow control simpler (no redundant code)
    def _writeout(do_we_slide):
        """
        Nested helper function that writes out the SNV count results for a window
        to the output file. The required arguments are a row object from the deque
        and a boolean value of whether or not to slide the window after writing out.
        """

        if len(window_deque) == 0:
            outbedg.write(_joinany("\t", [prev_chrom, win_start, win_end, 0]) + "\n")
        else:
            #count rare SNVs in window
            #snv_count = _count_rare_SNVs(window_deque, frequency)
            snv_count = _count_common_SNPs(window_deque, frequency)
            #write to file
            outbedg.write(_joinany("\t", [prev_chrom, win_start, win_end, snv_count]) + "\n")

        if do_we_slide == True:
            #sliding to next step
            _pop_bp_step(window_deque, win_step)



    def _append_OR_increment_and_append(row):
        """
        Nested helper function that either appends the current row to the deque
        or increments window to the row's position and appends to the deque, 
        writing out to the outputfile while doing so. The required argument is
        a row object from the deque.
        """

        #because these variables will be modified
        nonlocal win_start, win_end, win_step

        #decide whether to append to deque or increment window then append to deque
        if row.pos < win_start:
            exit("Record position before window start position")
        if row.pos <= win_end:
            window_deque.append(row)
        else:
            #increment window until position is reached, each window increment writes count to output file until current row position is reached
            for i in range(win_end, row.pos, window_step):

                #count SNVs and write results to file and slide the window by a step
                _writeout(True)

                #increment window
                win_end = i + window_step
                win_start = win_end - window_size + 1
                win_step = win_start + window_step

            #append current record to deque
            window_deque.append(row)





    #setting up deque for sliding window
    window_deque = deque()

    #open files
    with open(MAFcsv, "r") as csvf, open(f"{output}.bedgraph", "a") as outbedg:

        #write header line
        #outbedg.write("\t".join(["#CHROM", "START", "END", f"Num_{frequency}_SNVs"]) + "\n")
        outbedg.write("\t".join(["#CHROM", "START", "END", f"Num_{frequency}_SNPs"]) + "\n")

        #chromosome book keeping
        prev_chrom = None

        #window book keeping
        win_start = 1
        win_end = window_size
        win_step = window_step

        #loop through lines
        for line in csvf:
            #skip header line
            if line.startswith("CHROM"):
                continue

            #holding on current line
            row = MafCsvRowLoader(line)
            #setting the first chromosome, this only happens at the first row
            if prev_chrom is None:
                prev_chrom = row.chrom

            #check if row is on the same chromosome, if so you may proceed
            if row.chrom == prev_chrom:
                _append_OR_increment_and_append(row)
            else:
                #count SNVs, write out deque, and no need to slide window
                _writeout(False)

                #dump out deque
                window_deque.clear()

                #set new chromosome
                prev_chrom = row.chrom

                #window book keeping
                win_start = 1
                win_end = window_size
                win_step = window_step

                #decide whether to append to deque or increment window then append to deque
                _append_OR_increment_and_append(row)

        #clean up last window
        if len(window_deque) > 0:
            #count rare SNVs in window
            #snv_count = _count_rare_SNVs(window_deque, frequency)
            snv_count = _count_common_SNPs(window_deque, frequency)
            #write to file
            outbedg.write(_joinany("\t", [row.chrom, win_start, win_end, snv_count]) + "\n")

    print("Done!")






### RUNNING THE PROGRAM ###

# persite_MAF("beagle_biallelic_SNPs_1000GP30X_AllChr.vcf.gz")

# sliding_window_SNV_count("beagle_biallelic_SNPs_1000GP30X_AllChr_MAF.csv", 10000, 1000, 0.01, "win10Kb_step1Kb_0.01_SNV_count")
# sliding_window_SNV_count("beagle_biallelic_SNPs_1000GP30X_AllChr_MAF.csv", 100000, 10000, 0.01, "win100Kb_step10Kb_0.01_SNV_count")
# sliding_window_SNV_count("beagle_biallelic_SNPs_1000GP30X_AllChr_MAF.csv", 500000, 50000, 0.01, "win500Kb_step50Kb_0.01_SNV_count")


sliding_window_SNV_count("beagle_biallelic_SNPs_1000GP30X_AllChr_MAF.csv", 10000, 1000, 0.01, "win10Kb_step1Kb_0.01_SNP_count")
sliding_window_SNV_count("beagle_biallelic_SNPs_1000GP30X_AllChr_MAF.csv", 100000, 10000, 0.01, "win100Kb_step10Kb_0.01_SNP_count")
sliding_window_SNV_count("beagle_biallelic_SNPs_1000GP30X_AllChr_MAF.csv", 500000, 50000, 0.01, "win500Kb_step50Kb_0.01_SNP_count")