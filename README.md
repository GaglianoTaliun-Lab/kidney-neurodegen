# kidney-neurodegen

## Global genetic correlation analysis (LDSC)

1. GWAS summary statistics files found in the "/sumstats_data" directory are formatted in a format acceptable by LDSC (ie: have columns "CHR" and "BP to be able to annotate rsids) using scripts "formatting_GWASsumstats.R" and "VCF_to_sumstats_format.R" and saves them in the "/processed_sumstats" directory.
   
2. "Global_Genetic_Correlation_LDSC.sh" goes through all the different combinations using the new formatted GWAS summary statistics (/processed_sumstats directory) and perform munge_sumstats.py and ldsc.py and saves the results (/ldsc_results/ldscore_regression directory)
3. All the combination of phenotype's heritability is also found and saved (/ldsc_results/heritability) using "run_heritability.sh".

# LDSC result analysis
4. The script "heatmap_global_genetic_correlations.R" makes 3 heatmaps showing all the different global genetic correlations obtained. 

## Local genetic correlation analysis (LAVA)
# Prepare files for LAVA.

5. The script "input_info.R" uses a file "case_control_input_meta_analyzed.csv", previously made containing all the phenotype names, and number of cases and controls for all AD/PD, PD meta-analyzed (with and without proxies) in order to find the paths of each of those phenotypes and saves them in "input_info_PD_meta.txt".

6. The script "sample_overlap.R" links the results of cross-trait (LDSC) to LAVA.
   -> Step 1 : Uses all the results obtained from LDSC (in ldscore_regression directory) and combine them into 1 file and gives "ldsc_correlations_for_sample_overlap.txt" saved in       the directory "ldsc_corr".
   -> Step 2  : Creates the sample_overlap matrix for each correlation using "ldsc_correlations_for_sample_overlap.txt" (step 1 output) and saves them into the 
      "sample_overlap_files" directory.

# Run LAVA 
7. The script "LAVA.sh" goes through all "LAVA.R" for each combination of phenotypes to give univariate results with a Bonferroni thrsehold and Bivariate results.

# LAVA result analysis
8. The bivariate results obtained with LAVA are used by the script "unique_shared_loci.R" to find the number of shared and unique loci between proxy and non-proxy AD (sex-combined) and PD (sex-stratified, sex-combined and sex-combined meta-analyzed) for 3 different categories :
      Univariate test (passing the Bonferroni thrsehold) --> 
      Bivariate test (passing a p<0.05 thrsehold)
      Bivariate test (passing a Bonferroni thrsehold)
   
9. "venn_diagram_data.R" uses the results obtained by 'unique_shared_loci.R" (/LAVA/ven_diagram directory) to make the 
