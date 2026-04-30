#!/bin/bash


bedtools intersect -a hwe_0.01pvalue_windowed_density_100000bpwin_10000step.nocentromere.bedgraph -b vollger_hg38_SDs.merged.bed -wa > windows_intersecting_SDs.bed
bedtools intersect -a hwe_0.01pvalue_windowed_density_100000bpwin_10000step.nocentromere.bedgraph -b vollger_hg38_LargeSDs.merged.bed -wa > windows_intersecting_LargeSDs.bed
bedtools intersect -a hwe_0.01pvalue_windowed_density_100000bpwin_10000step.nocentromere.bedgraph -b vollger_hg38_SDs.merged.bed -v > windows_NOTintersecting_SDs.bed

bedtools intersect -a hwe_0.01pvalue_windowed_density_100000bpwin_10000step.nocentromere.bedgraph -b vollger_hg38_SDs.merged.bed -loj > windows_mapped_to_SDs.bed
bedtools intersect -a hwe_0.01pvalue_windowed_density_100000bpwin_10000step.nocentromere.bedgraph -b vollger_hg38_SDs.merged.bed -wao > windows_SD_overlaps.bed