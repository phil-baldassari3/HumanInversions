#!/bin/bash

set -e

###SUBSETTING###
echo "Subsetting VCF for selected peaks"

#chr1
bcftools view -r chr1:145906002-146806000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak1.chr1_145906002_146806000.invbrk1.chr1_146298110_NA.vcf.gz
bcftools view -r chr1:148183002-149682000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak2.chr1_148183002_149682000.invbrk1.chr1_NA_148672872.vcf.gz
bcftools view -r chr1:156040002-157280000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak3.chr1_156040002_157280000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr1:161030002-161990000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak4.chr1_161030002_161990000.invbrk0.NA_NA_NA.vcf.gz


#chr2
bcftools view -r chr2:76140002-78120000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak5.chr2_76140002_78120000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr2:86590002-88340000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak6.chr2_86590002_88340000.invbrk2.chr2_87987172_NA.vcf.gz
bcftools view -r chr2:111080002-112210000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak7.chr2_111080002_112210000.invbrk2.chr2_NA_111255403.vcf.gz
bcftools view -r chr2:120190002-121520000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak8.chr2_120190002_121520000.invbrk0.NA_NA_NA.vcf.gz


#chr3
bcftools view -r chr3:62950002-63710000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak9.chr3_62950002_63710000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr3:65180002-66220000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak10.chr3_65180002_66220000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr3:75080002-76170000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak11.chr3_75080002_76170000.invbrk0.NA_NA_NA.vcf.gz


#chr4
bcftools view -r chr4:14130002-15220000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak12.chr4_14130002_15220000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr4:25260002-26550000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak13.chr4_25260002_26550000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr4:39660002-41060000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak14.chr4_39660002_41060000.invbrk0.NA_NA_NA.vcf.gz


#chr5
bcftools view -r chr5:69190002-71960000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak15.chr5_69190002_71960000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr5:99160002-100350000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak16.chr5_99160002_100350000.invbrk0.NA_NA_NA.vcf.gz


#chr8
bcftools view -r chr8:3560002-4560000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak17.chr8_3560002_4560000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr8:6790002-8380000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak18.chr8_6790002_8380000.invbrk3.chr8_7301025_NA.vcf.gz
bcftools view -r chr8:11760002-12860000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak19.chr8_11760002_12860000.invbrk3.chr8_NA_12598379.vcf.gz
bcftools view -r chr8:16610002-17900000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak20.chr8_16610002_17900000.invbrk0.NA_NA_NA.vcf.gz


#chr10
bcftools view -r chr10:45963002-46972000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak21.chr10_45963002_46972000.invbrk4.chr10_46983452_NA.vcf.gz
bcftools view -r chr10:47340002-48190000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak22.chr10_47340002_48190000.invbrk4.chr10_NA_47468232.vcf.gz
bcftools view -r chr10:49770002-50350000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak23.chr10_49770002_50350000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr10:51700002-52560000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak24.chr10_51700002_52560000.invbrk0.NA_NA_NA.vcf.gz


#chr13
bcftools view -r chr13:35560002-36150000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak25.chr13_35560002_36150000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr13:36490002-36980000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak26.chr13_36490002_36980000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr13:45130002-47250000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak27.chr13_45130002_47250000.invbrk0.NA_NA_NA.vcf.gz


#chr15
bcftools view -r chr15:23010002-23361000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak28.chr15_23010002_23361000.invbrk5.chr15_23345460_NA.vcf.gz
bcftools view -r chr15:28226002-28784000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak29.chr15_28226002_28784000.invbrk5.chr15_NA_28389868.vcf.gz
bcftools view -r chr15:29990002-30910000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak30.chr15_29990002_30910000.invbrk6.chr15_30618104_NA.vcf.gz
bcftools view -r chr15:31950002-32820000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak31.chr15_31950002_32820000.invbrk6.chr15_NA_32153204.vcf.gz
bcftools view -r chr15:34250002-34800000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak32.chr15_34250002_34800000.invbrk0.NA_NA_NA.vcf.gz


#chr16
bcftools view -r chr16:16281002-16756000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak33.chr16_16281002_16756000.invbrk7.chr16_16721274_NA.vcf.gz
bcftools view -r chr16:18046002-18678000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak34.chr16_18046002_18678000.invbrk7.chr16_NA_18073542.vcf.gz
bcftools view -r chr16:20090002-20670000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak35.chr16_20090002_20670000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr16:21170002-22020000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak36.chr16_21170002_22020000.invbrk8.chr16_21583512_NA.vcf.gz
bcftools view -r chr16:22230002-22970000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak37.chr16_22230002_22970000.invbrk8.chr16_NA_22432381.vcf.gz


#chr17
bcftools view -r chr17:13680002-14470000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak38.chr17_13680002_14470000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr17:15310002-15950000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak39.chr17_15310002_15950000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr17:16380002-17090000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak40.chr17_16380002_17090000.invbrk9.chr17_16823491_NA.vcf.gz
bcftools view -r chr17:18270002-19250000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak41.chr17_18270002_19250000.invbrk9.chr17_NA_18384190.vcf.gz


#chr18
bcftools view -r chr18:9630002-9920000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak42.chr18_9630002_9920000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr18:9980002-10800000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak43.chr18_9980002_10800000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr18:11460002-12410000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak44.chr18_11460002_12410000.invbrk0.NA_NA_NA.vcf.gz


#chr19
bcftools view -r chr19:42650002-43820000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak45.chr19_42650002_43820000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr19:46200002-47180000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak46.chr19_46200002_47180000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr19:47610002-48180000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak47.chr19_47610002_48180000.invbrk10.chr19_47959661_NA.vcf.gz
bcftools view -r chr19:49910002-50280000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak48.chr19_49910002_50280000.invbrk10.chr19_NA_50091195.vcf.gz


#chr21
bcftools view -r chr21:26310002-26620000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak49.chr21_26310002_26620000.invbrk0.NA_NA_NA.vcf.gz
bcftools view -r chr21:27760002-28180000 /Users/philipbaldassari/Desktop/Phadnis_Lab/human_inversions/phased_vcfs/beagle_maf0.01_biallelic_SNPs_1000GP30X_AllChr.vcf.gz -Oz -o peak50.chr21_27760002_28180000.invbrk0.NA_NA_NA.vcf.gz



###INDEXING###
echo "Indexing output VCFs"


tabix -p vcf peak1.chr1_145906002_146806000.invbrk1.chr1_146298110_NA.vcf.gz
tabix -p vcf peak2.chr1_148183002_149682000.invbrk1.chr1_NA_148672872.vcf.gz
tabix -p vcf peak3.chr1_156040002_157280000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak4.chr1_161030002_161990000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak5.chr2_76140002_78120000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak6.chr2_86590002_88340000.invbrk2.chr2_87987172_NA.vcf.gz
tabix -p vcf peak7.chr2_111080002_112210000.invbrk2.chr2_NA_111255403.vcf.gz
tabix -p vcf peak8.chr2_120190002_121520000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak9.chr3_62950002_63710000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak10.chr3_65180002_66220000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak11.chr3_75080002_76170000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak12.chr4_14130002_15220000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak13.chr4_25260002_26550000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak14.chr4_39660002_41060000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak15.chr5_69190002_71960000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak16.chr5_99160002_100350000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak17.chr8_3560002_4560000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak18.chr8_6790002_8380000.invbrk3.chr8_7301025_NA.vcf.gz
tabix -p vcf peak19.chr8_11760002_12860000.invbrk3.chr8_NA_12598379.vcf.gz
tabix -p vcf peak20.chr8_16610002_17900000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak21.chr10_45963002_46972000.invbrk4.chr10_46983452_NA.vcf.gz
tabix -p vcf peak22.chr10_47340002_48190000.invbrk4.chr10_NA_47468232.vcf.gz
tabix -p vcf peak23.chr10_49770002_50350000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak24.chr10_51700002_52560000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak25.chr13_35560002_36150000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak26.chr13_36490002_36980000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak27.chr13_45130002_47250000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak28.chr15_23010002_23361000.invbrk5.chr15_23345460_NA.vcf.gz
tabix -p vcf peak29.chr15_28226002_28784000.invbrk5.chr15_NA_28389868.vcf.gz
tabix -p vcf peak30.chr15_29990002_30910000.invbrk6.chr15_30618104_NA.vcf.gz
tabix -p vcf peak31.chr15_31950002_32820000.invbrk6.chr15_NA_32153204.vcf.gz
tabix -p vcf peak32.chr15_34250002_34800000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak33.chr16_16281002_16756000.invbrk7.chr16_16721274_NA.vcf.gz
tabix -p vcf peak34.chr16_18046002_18678000.invbrk7.chr16_NA_18073542.vcf.gz
tabix -p vcf peak35.chr16_20090002_20670000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak36.chr16_21170002_22020000.invbrk8.chr16_21583512_NA.vcf.gz
tabix -p vcf peak37.chr16_22230002_22970000.invbrk8.chr16_NA_22432381.vcf.gz
tabix -p vcf peak38.chr17_13680002_14470000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak39.chr17_15310002_15950000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak40.chr17_16380002_17090000.invbrk9.chr17_16823491_NA.vcf.gz
tabix -p vcf peak41.chr17_18270002_19250000.invbrk9.chr17_NA_18384190.vcf.gz
tabix -p vcf peak42.chr18_9630002_9920000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak43.chr18_9980002_10800000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak44.chr18_11460002_12410000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak45.chr19_42650002_43820000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak46.chr19_46200002_47180000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak47.chr19_47610002_48180000.invbrk10.chr19_47959661_NA.vcf.gz
tabix -p vcf peak48.chr19_49910002_50280000.invbrk10.chr19_NA_50091195.vcf.gz
tabix -p vcf peak49.chr21_26310002_26620000.invbrk0.NA_NA_NA.vcf.gz
tabix -p vcf peak50.chr21_27760002_28180000.invbrk0.NA_NA_NA.vcf.gz

