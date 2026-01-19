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

- **UPGMA_and_hapcount_stats** Scripts used to test the statistical power of UPGMA and hapcount scans to identify inversions
    - `repeat_density_per_window.py`: Annotates each window with its corresponding repeat density (length of repeats / length of window) (repeats from UCSC RepeatMasker)
    - `branch_length_regressions.R`: Linear regression models of average branch length statistic with SNP density, repeat density, and presence of inversions
    - `hapcount_regressions.R`: Linear regression models of Haplotype Count statistic with SNP density and presence of inversions
    - `branch_len_permutation_testing.py`: Runs permutation Z-tests on Average Branch Length statistic and SNP density, repeat density, and presence of inversions
    - `hapcount_permutation_testing.py`: Runs permutation Z-tests on Haplotype Count statistic and SNP density and presence of inversions
    - `mark_inversion_windows.py`: annotates windows with labels denoting membership within an inversion feature (used for regression R scripts)
    - `mark_inversion_windows_with_bool.py`: annotates windows with 0 or 1 denoting membership within an inversion feature (used for ROC curves)
    - `inversion_classification_ROC_curves.py`: ROC curves showing the power of each statistic in classifying windows as inversions or not
    