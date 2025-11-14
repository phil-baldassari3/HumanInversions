#!/bin/bash/

# #inversion regions
# bcftools view -r chr9:112105263-114119310 -Oz -o chr9_1MB_inv_1MB_2014048bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz
# bcftools view -r chr19:20647331-23062459 -Oz -o chr19_1MB_inv_1MB_2415129bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz
# bcftools view -r chr17:35357258-38957978 -Oz -o chr17_1MB_inv_1MB_3600721bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz
# bcftools view -r chr15:29077909-33607507 -Oz -o chr15_1MB_inv_1MB_4529599bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz
# bcftools view -r chr15:21770522-29852548 -Oz -o ch15_1MB_inv_1MB_8082027bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz

# #control regions
# bcftools view -r chr9:81439687-83453735 -Oz -o chr9_control1_2014048bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz
# bcftools view -r chr9:1410271-3424319 -Oz -o chr9_control2_2014048bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr9.vcf.gz

# bcftools view -r chr19:11681438-14096567 -Oz -o chr19_control1_2415129bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz
# bcftools view -r chr19:54952346-57367475 -Oz -o chr19_control2_2415129bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz

# bcftools view -r chr17:9256449-12857170 -Oz -o chr17_control1_3600721bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz
# bcftools view -r chr17:50049695-53650416 -Oz -o chr17_control2_3600721bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr17.vcf.gz

# bcftools view -r chr15:66102101-70631700 -Oz -o chr15_control1_4529599bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz
# bcftools view -r chr15:49153208-53682807 -Oz -o chr15_control2_4529599bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz

# bcftools view -r chr15:77993044-86075071 -Oz -o chr15_control1_8082027bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz
# bcftools view -r chr15:35391135-43473162 -Oz -o chr15_control2_8082027bp.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr15.vcf.gz




bcftools view -r chr19:700000-700100 -Oz -o 100bp_chunk.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz
bcftools view -r chr19:700000-701000 -Oz -o 1000bp_chunk.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz
bcftools view -r chr19:700000-710000 -Oz -o 10000bp_chunk.vcf.gz beagle_phased_biallelic_SNPs_1000GP30X_chr19.vcf.gz