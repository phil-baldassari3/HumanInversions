#!/bin/bash


bcftools annotate --set-id '%POS' beagle_phased_biallelic_SNPs_1000GP30X_Autosomes.vcf.gz -Oz -o beagle_phased_biallelic_SNPs_1000GP30X_Autosomes.ann.vcf.gz