#!/bin/bash

java -Xmx20g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr15.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr15 map=~/bioinformatics_tools/beagle5/new_chr15.map nthreads=8
java -Xmx20g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr14.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr14 map=~/bioinformatics_tools/beagle5/new_chr14.map nthreads=8
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr6.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr6 map=~/bioinformatics_tools/beagle5/new_chr6.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr4.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr4 map=~/bioinformatics_tools/beagle5/new_chr4.map nthreads=10
java -Xmx64g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr1.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr1 map=~/bioinformatics_tools/beagle5/new_chr1.map nthreads=10


java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr19.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr19 map=~/bioinformatics_tools/beagle5/new_chr19.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr17.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr17 map=~/bioinformatics_tools/beagle5/new_chr17.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr16.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr16 map=~/bioinformatics_tools/beagle5/new_chr16.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr10.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr10 map=~/bioinformatics_tools/beagle5/new_chr10.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr9.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr9 map=~/bioinformatics_tools/beagle5/new_chr9.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr8.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr8 map=~/bioinformatics_tools/beagle5/new_chr8.map nthreads=10

java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr5.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr5 map=~/bioinformatics_tools/beagle5/new_chr5.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr11.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr11 map=~/bioinformatics_tools/beagle5/new_chr11.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr13.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr13 map=~/bioinformatics_tools/beagle5/new_chr13.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chrX.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chrX map=~/bioinformatics_tools/beagle5/new_chrX.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chrY.vcf.gz out=beagle_imputed_biallelic_SNPs_1000GP30X_chrY map=~/bioinformatics_tools/beagle5/new_chrY.map nthreads=10



java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr7.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr7 map=~/bioinformatics_tools/beagle5/new_chr7.map nthreads=10
java -Xmx64g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr3.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr3 map=~/bioinformatics_tools/beagle5/new_chr3.map nthreads=10


java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chrY.vcf.gz out=beagle_imputed_biallelic_SNPs_1000GP30X_chrY nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_nonPAR_chrX.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_nonPAR_chrX map=~/bioinformatics_tools/beagle5/new_chrX.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_PAR_chrX.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_PAR_chrX map=~/bioinformatics_tools/beagle5/new_chrX.map nthreads=10
java -Xmx64g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr2.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr2 map=~/bioinformatics_tools/beagle5/new_chr2.map nthreads=10




java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr12.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr12 map=~/bioinformatics_tools/beagle5/new_chr12.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr18.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr18 map=~/bioinformatics_tools/beagle5/new_chr18.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr20.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr20 map=~/bioinformatics_tools/beagle5/new_chr20.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr21.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr21 map=~/bioinformatics_tools/beagle5/new_chr21.map nthreads=10
java -Xmx32g -jar ~/bioinformatics_tools/beagle5/beagle.27Feb25.75f.jar gt=vcfs_1000GP_30X/full_biallelic_vcfs/biallelic_SNPs_1000GP30X_chr22.vcf.gz out=beagle_phased_biallelic_SNPs_1000GP30X_chr22 map=~/bioinformatics_tools/beagle5/new_chr22.map nthreads=10
