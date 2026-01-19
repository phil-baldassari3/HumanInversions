import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc
import matplotlib.pyplot as plt


def plot_ROC(df, parameter_column, truth_column):
    """
    Plots a ROC curve to assess the accuracy of using a given parameter to predict the identity of a genomic window (binary choice)
    
    :param df: windowed average branch length and other data parameters and groundtruth columns denoting whether a window is in or out (1 or 0) of a particular feature
    :param parameter_column: column name of the parameter you are assessing as a classifier (simple threshold model)
    :param truth_column: column name of 1s and 0s that show whether the window in in or out of that particular feature (groundtruth)
    """

    #classifier parameter
    vals = df[parameter_column].values

    #groundtruth classification
    truth = df[truth_column].values

    #run ROC
    fpr, tpr, thresholds = roc_curve(truth, vals)
    roc_auc = auc(fpr, tpr)

    #plotting
    title = f"ROC of {parameter_column} on {truth_column}"
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

    plt.savefig(f"{title}.png")
    plt.clf()


#hapcount_SNPwindow/cleaned_autosomes_SNPwindow1000_SNPstep500_hapcount_snpden_labeled_with_bool
data = pd.read_csv("hapcount_BPwindow/cleaned_autosomes_BPwindow3000_BPstep1500_hapcount_snpden_labeled_with_bool.tsv", sep="\t")

plot_ROC(data, "NUM_uniq_haps", "inversion")
plot_ROC(data, "NUM_uniq_haps", "large_inversion")
plot_ROC(data, "NUM_uniq_haps", "inv_breakpoint")
plot_ROC(data, "NUM_uniq_haps", "large_inv_breakpoint")