import random


#chromosome map: the rough beginning and end coordinates for variants on this chromosome
chr_map = {
    "chr9":(10123,138334684),
    "chr19":(60377,58607603),
    "chr17":(60104,83247307),
    "chr15":(17000095,101981147),
}

#inversion map: used to check whether generated region overlaps with an inversion region
inv_map = {
    "chr9":[(112105263, 114119310)],
    "chr19":[(20647331, 23062459)],
    "chr17":[(35357258, 38957978)],
    "chr15":[(29077909, 33607507),(21770522, 29852548)],
}


#targets: list of tuples of the form (chrom, length) used to generate randomized region coordinates
targets = [
    ("chr9", 2014048),
    ("chr19", 2415129),
    ("chr17", 3600721),
    ("chr15", 4529599),
    ("chr15", 8082027),
]

#how many randomized regions do you want for each target
rand_per_target = 2

#set seed
random.seed(563)



for t in targets:
    chrom = t[0]
    target_len = t[1]
    possible_start = chr_map[t[0]][0]
    possible_end = chr_map[t[0]][1] - target_len
    inv_to_avoid = inv_map[chrom]

    for i in range(rand_per_target):

        rand_start_coord = random.randint(possible_start, possible_end)
        rand_end_coord = rand_start_coord + target_len
        flag = ""

        #check if output coordinates are within the inversion region (will be flagged in output)
        for inv in inv_map[chrom]:
            if (rand_end_coord < inv[0]) or (rand_start_coord > inv[1]):
                continue
            else:
                flag = "OUTPUT WITHIN INVERSION REGION!"

        #print output to stdout
        print(f"{chrom}:{str(rand_start_coord)}-{str(rand_end_coord)}\t{flag}")



