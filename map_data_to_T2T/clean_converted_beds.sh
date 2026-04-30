#!/bin/bash

bedtools sort -i TEMP_BPwindow100000_BPstep10000_hapcount_converted_to_T2T.bedgraph | cut -f 1,2,3,4 > BPwindow100000_BPstep10000_hapcount_converted_to_T2T.bedgraph

bedtools sort -i TEMP_vollger_hg38_SDs.merged_converted_to_T2T.bed | cut -f 1,2,3 > vollger_hg38_SDs.merged_converted_to_T2T.bed
bedtools sort -i TEMP_lowessfindpeaks90_BPwindow100000_BPstep10000_hapcount.merged_converted_to_T2T.bed | cut -f 1,2,3 > lowessfindpeaks90_BPwindow100000_BPstep10000_hapcount.merged_converted_to_T2T.bed
bedtools sort -i TEMP_peaks_0inv_0SDs_0SDbp_converted_to_T2T.bed | cut -f 1,2,3 > peaks_0inv_0SDs_0SDbp_converted_to_T2T.bed
bedtools sort -i TEMP_peaks_1inv_5SDs_5000SDbp_converted_to_T2T.bed | cut -f 1,2,3 > peaks_1inv_5SDs_5000SDbp_converted_to_T2T.bed
