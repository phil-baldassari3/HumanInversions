# HumanInversions

Scripts for project to identify genomic inversions using human genotype data.

Directories and Scripts:

- **vcf_processing_and_phasing**: Scripts used to:
    1. download 30X genotype vcfs from 1000 Genomes Project
    2. filter for biallelic SNVs only and index with tabix
    3. phase and impute missing genotypes using Beagle5.5
    - note that the PAR and nonPAR regions on X were phased spearately
    - note that the Y chromosome was imputed using Beagle and bcftools was used to remove female samples

- **benchmarking_beagle**: Scripts used benchmark beagle on my system

- **het_and_maf_scans**: Scripts used to run genome scans of Heterozygosity and Minor Allele Frequency
    - These scans were not able to identify inversions

- **PCA_scan**: Scripts used to run a windowed Principle Component scan on haplotypes across the genome
    - These scans were not able to identify inversions

- **first_round_UPGMA_scan**: Scripts used to test feasibility of running a windowed UPGMA scan on haplotypes
    - This method showed promise in identifying inversions breakpoints

- **UPGMA_scan**: Scripts to run a windowed UPGMA scan
    - `SNP_windows_UPGMA_windowed_scan.py`: runs UPGMA in defined SNP-windows and SNP-steps
    - `bp_windows_UPGMA_windowed_scan.py`: runs UPGMA in defined bp-windows and bp-steps

- **hapcount_scan**: Scripts to run a scan of counting the number of unique haplotypes in each window (faster than UPGMA)
    - `SNPwindow_hap_counter.py`: counts haplotypes in defined SNP-windows and SNP-steps
    - `BPwindow_hap_counter.py`: counts haplotypes in defined SNP-windows and bp-steps
    - `hapcount_scan_v1.py`: first version of haplotype count scan in SNP windows

- **UPGMA_and_hapcount_stats**: Scripts used to test the statistical power of UPGMA and hapcount scans to identify inversions
    - `repeat_density_per_window.py`: Annotates each window with its corresponding repeat density (length of repeats / length of window) (repeats from UCSC RepeatMasker)
    - `branch_length_regressions.R`: Linear regression models of average branch length statistic with SNP density, repeat density, and presence of inversions
    - `hapcount_regressions.R`: Linear regression models of Haplotype Count statistic with SNP density and presence of inversions
    - `branch_len_permutation_testing.py`: Runs permutation Z-tests on Average Branch Length statistic and SNP density, repeat density, and presence of inversions
    - `hapcount_permutation_testing.py`: Runs permutation Z-tests on Haplotype Count statistic and SNP density and presence of inversions
    - `mark_inversion_windows.py`: annotates windows with labels denoting membership within an inversion feature (used for regression R scripts)
    - `mark_inversion_windows_with_bool.py`: annotates windows with 0 or 1 denoting membership within an inversion feature (used for ROC curves)
    - `inversion_classification_ROC_curves.py`: ROC curves showing the power of each statistic in classifying windows as inversions or not

- **variant_density_scan**: Scripts to perform a genome scan of a count of all variants per window or variants above or below 0.01 frequency
    - `low_freq_SNV_scan.py`: Counts number of SNVs with allele frequency less than 0.01 in windows
    - `windowed_variance_scan.py`: Counts the variance in SNV counts across multiple windows
    - `absZscore_transform.py`: Absolute value Z-Score transform of variant counts per window
    - `compute_all_bialleleic_variant_count.py`: Counts number of SNVs in windows
    - `compute_SNP_to_SNV_ratio.py`: Computes the ratio of SNVs with greater than 0.01 frequency to SNVs with less than 0.01 frequency in windows

- **testing_hapcount**: Scripts to test the ability of windowed Haplotype Proportion to predict Inversions and Segmental Duplications
    - `clean_hapcount_bedgraph.py`: Clips out centromeres from Haplotype Proportion data
    - `mark_InvBrk_SD_windows_with_bool.py`: marks windows that map to inversion breakpoints or SDs with a 0 or 1
    - `mark_InvBrk_SD_windows.py`: marks windows that map to inversion breakpoints or SDs with a label
    - `hapcount_permutation_testing.py`: Permutation test between Haplotype Proportion of Inv/SD windows vs. whole genomes
    - `hapcount_regressions.R`: Haplotype Proportion regressions and other stats
    - `inversion_classification_ROC_curves.py`: ROC curve of Haplotype Proportion as a predictor of Invs/SDs
    - `peak_finding.py`: Finds high Haplotype Proportion peaks
    - `counting_intersects.sh`: Runs bedtools to find intersects between called peaks and Invs/SDs
    - `plot_mapping_proportions.R`: Plots how many peaks mapped to Invs/SDs

- **long_range_LD**: Scripts to estimate bootstrapped Linkage-Disequilibrium betwen called peaks
    - `bootstrap_LD_v1.1.py`: version 1 of script that runs LD using plink between two regions, this version does not bootstrap
    - `bootstrap_LD.py`: runs bootstrap LD using plink between pairs of regions
    - `dprime_gmeans_plotting.R`: plotting non-bootstrapped LD analysis
    - `individual_LD_stats.py`: computes "LDstat" for each individual per comparision
    - `plot_bootstrapLD.R`: plotting bootstrapped LD analysis
    - `plot_pairwise_region_LD.R`: plotting D-prime between two regions
    - `subset_vcf_for_LD_testing.sh`: subset vcfs for comparisions

- **inverted_homology_qunatification**: Scripts that test whether the amopunt of inverted homology between regions can be used to predict inversions
    - `invdup_permutation_test.py`: permutation test of the amount of inverted homology around inversions vs. whole genome
    - `ROC_invdup_predict_invs.py`: ROC curve of amoutn of inverted homology between regions to predict inversions
    - `find_windowed_SD_density.py`: finds windowed density of SDs from the hg38 reference SDs callset from Vollger et al. 2022
    - `find_SD_density_peaks.py`: finds peaks of high (greater than 25%) SD density

- **hwe_scan**: Scripts to see if high Haplotype Proportion is just a result of bad read mapping/SNPs out of HWE
    - `label_SNPs_with_POS.sh`: adds a SNPid column to a vcf that just labels each SNP with its genomic position (needed because plink does not output position)
    - `run_plink_hwe.sh`: runs plink to test for HWE for each variant
    - `hwe_windowed_p_density.py`: Scans for the density of significant p-values in sliding windows
    - `intersects_of_windows_and_SDs.sh`: maps windowed sigP-value data to SDs
    - `HWE_pval_density_and_SDs.R`: plots results of sigP-value density and runs stats

- **map_data_to_T2T**: Scripts that map data to T2T reference
    - `download_HGSVC3_SDs.sh`: Downloads SD data from HGSVC3
    - `clean_dedup_SD_bed.py`: sorts and dedupes downloaded SDs
    - `clean_converted_beds.sh`: selects columns to make proper beds and bedgraphs
    - `cat_sort_uniq_sample_SDs.sh`: merges SDs from population data and sorts and dedups


