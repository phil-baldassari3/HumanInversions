import pandas as pd
import matplotlib.pyplot as plt


inv_map = {
    "chr1":[(205209479, 205209706), (187495697, 187497600)],
    "chr4":[(40233407, 40235440), (3894504, 9763113)],
    "chr6":[(130527041, 130531150), (31041443, 31042319)],
    "chr8":[(7064966, 12716088)],
    "chr9":[(113105263, 113119310)],
    "chr10":[(36819152, 58318428)],
    "chr14":[(60604530, 60613248), (34540845, 34562270)],
    "chr15":[(22770522, 28852548), (30077909, 32607507)],
    "chr17":[(36357258, 37957978)],
    "chr19":[(21647331, 22062459)]
    }


def hetplot(csvfile, chrom, window):

    df = pd.read_csv(csvfile)

    #groundtruth inversions
    grndt_invs = inv_map[chrom]

    plt.figure(figsize=(15, 5))
    if window >= 1000:
        for i in df.columns[2:-2]:
            plt.scatter(df["POS"], df[i], s=0.5, color="gray", alpha=0.2)
    plt.plot(df["POS"], df["p10_mean"], color='blue', label="p^10 Generalized Mean")
    plt.plot(df["POS"], df["Geometric_mean"], color="orange", label="Geometric Mean")
    #plt.yscale('log')
    plt.xlabel("chr15")
    plt.ylabel("Heterozygosity")
    plt.ticklabel_format(style='plain')

    for span in grndt_invs:
        plt.axvspan(span[0], span[1], color='red', alpha=0.5)

    plt.legend()
    plt.title(f"{chrom} Heterozygosity {window} SNP windows")
    plt.savefig(f"Het_{chrom}_{window}_SNP_windows.png")
    #plt.show()
    plt.clf()




hetplot("Heterozygosity_chr19_100_SNP_windows.csv", "chr19", 100)
hetplot("Heterozygosity_chr19_1000_SNP_windows.csv", "chr19", 1000)
hetplot("Heterozygosity_chr19_10000_SNP_windows.csv", "chr19", 10000)

hetplot("Heterozygosity_chr17_100_SNP_windows.csv", "chr17", 100)
hetplot("Heterozygosity_chr17_1000_SNP_windows.csv", "chr17", 1000)
hetplot("Heterozygosity_chr17_10000_SNP_windows.csv", "chr17", 10000)

hetplot("Heterozygosity_chr15_100_SNP_windows.csv", "chr15", 100)
hetplot("Heterozygosity_chr15_1000_SNP_windows.csv", "chr15", 1000)
hetplot("Heterozygosity_chr15_10000_SNP_windows.csv", "chr15", 10000)

hetplot("Heterozygosity_chr14_100_SNP_windows.csv", "chr14", 100)
hetplot("Heterozygosity_chr14_1000_SNP_windows.csv", "chr14", 1000)
hetplot("Heterozygosity_chr14_10000_SNP_windows.csv", "chr14", 10000)