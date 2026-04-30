#!/bin/bash

set -e

# #SDs intersecting Invs
# bedtools intersect -a inv_SD_bed_files/vollger_hg38_SDs.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > SDs_count_mapping_Invs.bed

# #Invs intersecting SDs
# bedtools intersect -a inv_SD_bed_files/porubsky_inversions.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Invs_count_mapping_SDs.bed

#SDs intersecting InvBrks
bedtools intersect -a inv_SD_bed_files/vollger_hg38_SDs.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > SDs_count_mapping_InvBrks.bed

#InvBrks intersecting SDs
bedtools intersect -a inv_SD_bed_files/porubsky_inversions_breakpoints.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > InvBrks_count_mapping_SDs.bed


# ###Peak finding methods intersecting SDs
# bedtools intersect -a putative_peaks/lowess90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_lowess90_100Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/lowess95_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_lowess95_100Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/findpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_findpeaks90_100Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/findpeaks95_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_findpeaks95_100Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_lowessfindpeaks90_100Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/lowess90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_lowess90_10Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/lowess95_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_lowess95_10Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/findpeaks90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_findpeaks90_10Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/findpeaks95_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_findpeaks95_10Kbwin_count_mapping_SDs.bed
# bedtools intersect -a putative_peaks/lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/vollger_hg38_SDs.merged.bed -c > Peaks_lowessfindpeaks90_10Kbwin_count_mapping_SDs.bed



# ###Peak finding methods intersecting Invs
# bedtools intersect -a putative_peaks/lowess90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_lowess90_100Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/lowess95_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_lowess95_100Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/findpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_findpeaks90_100Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/findpeaks95_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_findpeaks95_100Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_lowessfindpeaks90_100Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/lowess90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_lowess90_10Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/lowess95_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_lowess95_10Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/findpeaks90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_findpeaks90_10Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/findpeaks95_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_findpeaks95_10Kbwin_count_mapping_Invs.bed
# bedtools intersect -a putative_peaks/lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions.bed -c > Peaks_lowessfindpeaks90_10Kbwin_count_mapping_Invs.bed




###Peak finding methods intersecting InvBrks
bedtools intersect -a putative_peaks/lowess90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_lowess90_100Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/lowess95_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_lowess95_100Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/findpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_findpeaks90_100Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/findpeaks95_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_findpeaks95_100Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow100000_BPstep10000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_lowessfindpeaks90_100Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/lowess90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_lowess90_10Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/lowess95_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_lowess95_10Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/findpeaks90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_findpeaks90_10Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/findpeaks95_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_findpeaks95_10Kbwin_count_mapping_InvBrks.bed
bedtools intersect -a putative_peaks/lowessfindpeaks90_Cen_Clipped_ALL_CHROM_BPwindow10000_BPstep1000_hapcount.merged.bed -b inv_SD_bed_files/porubsky_inversions_breakpoints.bed -c > Peaks_lowessfindpeaks90_10Kbwin_count_mapping_InvBrks.bed


