import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt


def plot_ROC(df, parameter_column, truth_column, filelab):
    """
    Plots a ROC curve to assess the accuracy of using a given parameter to predict the identity of a genomic window (binary choice)
    
    :param df: windowed average branch length and other data parameters and groundtruth columns denoting whether a window is in or out (1 or 0) of a particular feature
    :param parameter_column: column name of the parameter you are assessing as a classifier (simple threshold model)
    :param truth_column: column name of 1s and 0s that show whether the window in in or out of that particular feature (groundtruth)
    :param filelab: label used for saving a unique filename e.g. "10Kb_windows"
    """

    #classifier parameter
    vals = df[parameter_column].values

    #groundtruth classification
    truth = df[truth_column].values

    #run ROC
    fpr, tpr, thresholds = roc_curve(truth, vals)
    roc_auc = auc(fpr, tpr)

    #plotting
    title = f"ROC of {parameter_column} on {truth_column}".replace("Proportion_of_Unique_Haplotypes", "Unique Haplotype Proportion")
    plt.figure(figsize=(6,6))
    plt.plot(fpr, tpr, label=f"AUC = {round(roc_auc, 3)}")
    plt.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.7, label="No Skill")
    plt.axvline(x=0.05, color="black", linewidth=0.7, label="5% FPR cutoff")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    #plt.show()

    plt.savefig(f"{title}_{filelab}.png")
    plt.clf()



def _run_ROC(values, truths):
    """
    Helper function that computes ROC on an array of values adn truth values and
    returns a tuple with a fpr array, tpr array, and auc
    """

    #run ROC
    fpr, tpr, thresholds = roc_curve(truths, values)
    roc_auc = auc(fpr, tpr)

    return fpr, tpr, roc_auc


def plot_ROCv2(df, parameter_column, truth_columns, filelab):
    """
    Plots a ROC curve to assess the accuracy of using a given parameter to predict the identity of a genomic window (binary choice)
    
    :param df: windowed average branch length and other data parameters and groundtruth columns denoting whether a window is in or out (1 or 0) of a particular feature
    :param parameter_column: column name of the parameter you are assessing as a classifier (simple threshold model)
    :param truth_column: column name of 1s and 0s that show whether the window in in or out of that particular feature (groundtruth)
    :param filelab: label used for saving a unique filename e.g. "10Kb_windows"
    """

    #classifier parameter
    vals = df[parameter_column].values

    #plotting
    title = f"ROC of {parameter_column}".replace("Proportion_of_Unique_Haplotypes", "Unique Haplotype Proportion")
    plt.figure(figsize=(6,6))

    for col in truth_columns:

        #groundtruth classification
        truth = df[col].values
        outfpr, outtpr, out_auc = _run_ROC(vals, truth)
        plt.plot(outfpr, outtpr, label=f"{col} (AUC = {round(out_auc, 3)})")


    plt.plot([0, 1], [0, 1], linestyle="--", color="black", linewidth=0.7, label="No Skill")
    plt.axvline(x=0.05, color="black", linewidth=0.7, label="5% FPR cutoff")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    #plt.show()

    plt.savefig(f"{title}_{filelab}.png")
    plt.clf()




def FDR_vs_TPR(df, parameter_column, truth_columns, filelab):
    """
    Function plots the TPR vs. FDR of the unique haplotype proportion to predict inversions and SDs.

    :param df: windowed average branch length and other data parameters and groundtruth columns denoting whether a window is in or out (1 or 0) of a particular feature
    :param parameter_column: column name of the parameter you are assessing as a classifier (simple threshold model)
    :param truth_column: column name of 1s and 0s that show whether the window in in or out of that particular feature (groundtruth)
    :param filelab: label used for saving a unique filename e.g. "10Kb_windows"
    """

    #plot
    plt.figure(figsize=(6,6))

    #get list of thresholds
    threshs = list(set(df[parameter_column].to_list()))
    threshs.sort(reverse=True)

    #iterate through thruth columns
    for col in truth_columns:

        total_tp = df[col].sum()
        tpr = []
        fdr = []

        #iterate through thresholds
        for thresh in threshs:

            filtered_df = df[df[parameter_column] >= thresh]
        
            tp = filtered_df[col].sum()
            fp = len(filtered_df) - filtered_df[col].sum()
            fp_tp = len(filtered_df)

            tpr.append(tp / total_tp)
            fdr.append(fp / fp_tp)


        plt.plot(tpr, fdr, label=f"{col}")

    plt.ylabel("False Discovery Rate")
    plt.xlabel("True Positive Rate")
    plt.title("False Discovery Rate vs. True Positive Rate\nof Unique Haplotype Proportion Threshold")
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(f"{filelab}_False Discovery Rate vs. True Positive Rate of Unique Haplotype Proportion Threshold.png")
    plt.clf()





def FDR_vs_Threshold(df, parameter_column, truth_columns, filelab):
    """
    Function plots the FFR vs. Threshold of the unique haplotype proportion to predict inversions and SDs.

    :param df: windowed average branch length and other data parameters and groundtruth columns denoting whether a window is in or out (1 or 0) of a particular feature
    :param parameter_column: column name of the parameter you are assessing as a classifier (simple threshold model)
    :param truth_column: column name of 1s and 0s that show whether the window in in or out of that particular feature (groundtruth)
    :param filelab: label used for saving a unique filename e.g. "10Kb_windows"
    """

    #plot
    plt.figure(figsize=(8,6))

    #get list of thresholds
    threshs = list(set(df[parameter_column].to_list()))
    threshs.sort(reverse=True)

    #iterate through thruth columns
    for col in truth_columns:

        total_tp = df[col].sum()
        tpr = []
        fdr = []

        #iterate through thresholds
        for thresh in threshs:

            filtered_df = df[df[parameter_column] >= thresh]
        
            tp = filtered_df[col].sum()
            fp = len(filtered_df) - filtered_df[col].sum()
            fp_tp = len(filtered_df)

            tpr.append(tp / total_tp)
            fdr.append(fp / fp_tp)

        plt.plot(threshs, fdr, label=f"{col}")

    plt.ticklabel_format(style='plain', axis='x')
    plt.xlabel("Threshold (Proportion of Unique Haplotypes)")
    plt.ylabel("False Discovery Rate")
    plt.title("False Discovery Rate vs. Proportion of Unique Haplotypes as a Threshold")
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(f"{filelab}_False Discovery Rate vs. Proportion of Unique Haplotypes as a Threshold.png")
    plt.clf()



# ###100Kb
# data = pd.read_csv("Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount_boollabeled.tsv", sep="\t")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "InvBrk", "100Kb")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "LargeInvBrk", "100Kb")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "SD", "100Kb")


# ###10Kb
# data = pd.read_csv("Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount_boollabeled.tsv", sep="\t")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "InvBrk", "10Kb")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "LargeInvBrk", "10Kb")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "SD", "10Kb")


# ###10KSNP
# data = pd.read_csv("Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount_boollabeled.tsv", sep="\t")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "InvBrk", "10KSNP")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "LargeInvBrk", "10KSNP")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "SD", "10KSNP")


# ###1KSNP
# data = pd.read_csv("Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount_boollabeled.tsv", sep="\t")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "InvBrk", "1KSNP")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "LargeInvBrk", "1KSNP")
# plot_ROC(data, "Proportion_of_Unique_Haplotypes", "SD", "1KSNP")







###100Kb
data = pd.read_csv("Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount_boollabeled.tsv", sep="\t")
#plot_ROCv2(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "100Kb")
FDR_vs_TPR(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "100Kb")
FDR_vs_Threshold(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "100Kb")


###10Kb
#data = pd.read_csv("Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount_boollabeled.tsv", sep="\t")
#plot_ROCv2(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "10Kb")




###10KSNP
data = pd.read_csv("Cen_Clipped_All_CHROM_SNPwindow10000_SNPstep1000_hapcount_boollabeled.tsv", sep="\t")
#plot_ROCv2(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "10KSNP")
FDR_vs_TPR(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "10KSNP")
FDR_vs_Threshold(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "10KSNP")




###1KSNP
#data = pd.read_csv("Cen_Clipped_All_CHROM_SNPwindow1000_SNPstep100_hapcount_boollabeled.tsv", sep="\t")
#plot_ROCv2(data, "Proportion_of_Unique_Haplotypes", ["InvBrk", "LargeInvBrk", "SD", "LargeSD"], "1KSNP")

