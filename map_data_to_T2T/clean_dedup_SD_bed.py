import os


def clean_dedup_SD_bedfile(bedfile):
    """
    Function takes in a bedfile output from SEDEF with called Segmental Duplications
    and saves a new bedfile that lists each unique SD from the origianl output columns
    1-3 (CHROM, START, END)

    :param bedfile: filename of SEDEF output
    :type bedfile: str
    """

    print(f"Cleaning {bedfile}...")

    #setting empty list to load with SD spans
    sds = []

    #open input file and iterate through lines
    with open(bedfile, "r") as inputbed:
        for line in inputbed:

            #skip header
            if line.startswith("#"):
                continue

            #grab SD span
            linelist = line.strip().split("\t")
            chrom = linelist[0].replace("_RagTag", "")
            #only keep SDs that mapped to chromosomes
            if "chr" not in chrom:
                continue
            start = int(linelist[1])
            end = int(linelist[2])

            #append to list
            sds.append((chrom, start, end))

    #sort and dedup list of sds
    unique_sds = sorted(list(set(sds)))

    #write to output
    with open(bedfile.replace(".bed", ".CleanDedup.bed"), "w") as outputbed:
        for c, s, e in unique_sds:
            outputbed.write(f"{c}\t{s}\t{e}\n")

    print(f"{bedfile} cleaning complete!\n")




#Run on all bedfiles in directory
for file in os.listdir():
    if file.endswith(".bed"):
        clean_dedup_SD_bedfile(file)



