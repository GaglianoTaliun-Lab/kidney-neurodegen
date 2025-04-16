# kidney-neurodegen

This repository hosts the code for the work described in:

Sadaf Gawhary, Le Chang, Lyza Maameri, Wiame Belbellaj, Frida Lona-Durazo, Sarah A Gagliano Taliun (2025) Assessment of global and local genetic correlations between kidney-related traits and late-onset neurodegenerative diseases.

## Global genetic correlation analysis (LDSC)

1. GWAS summary statistics files found in the "/sumstats_data" directory are formatted in a format acceptable by LDSC (ie: have columns "CHR" and "BP to be able to annotate rsids) using scripts "formatting_GWASsumstats.R" and "VCF_to_sumstats_format.R" and saves them in the "/processed_sumstats" directory.
   
2. "Global_Genetic_Correlation_LDSC.sh" goes through all the different combinations using the new formatted GWAS summary statistics (/processed_sumstats directory) and perform munge_sumstats.py and ldsc.py and saves the results (/ldsc_results/ldscore_regression directory)
   
3. All the combination of phenotype's heritability is also found and saved (/ldsc_results/heritability) using "run_heritability.sh".


LDSC result analysis
4. The script "heatmap_global_genetic_correlations.R" makes 3 heatmaps showing all the different global genetic correlations obtained. 

## Local genetic correlation analysis (LAVA)
Prepare files for LAVA.

5. The script "input_info.R" uses a file "case_control_input_meta_analyzed.csv", previously made containing all the phenotype names, and number of cases and controls for all AD/PD, PD meta-analyzed (with and without proxies) in order to find the paths of each of those phenotypes and saves them in "input_info_PD_meta.txt".
   
6. The script "sample_overlap.R" links the results of cross-trait (LDSC) to LAVA.
   
   -> Step 1 : Uses all the results obtained from LDSC (in ldscore_regression directory) and combine them into 1 file and gives "ldsc_correlations_for_sample_overlap.txt" saved in       the directory "ldsc_corr".
   
   -> Step 2  : Creates the sample_overlap matrix for each correlation using "ldsc_correlations_for_sample_overlap.txt" (step 1 output) and saves them into the 
      "sample_overlap_files" directory.


Run LAVA 

7. The script "LAVA.sh" goes through all "LAVA.R" for each combination of phenotypes to give univariate results with a Bonferroni thrsehold and Bivariate results.

LAVA result analysis

8. The bivariate results obtained with LAVA are used by the script "venn_diagram_data.R" to find the number of shared and unique loci between proxy and non-proxy AD (sex-combined) and PD (sex-stratified, sex-combined and sex-combined meta-analyzed) for 3 different categories :
      Univariate test (passing the Bonferroni thrsehold)
      Bivariate test (passing a p<0.05 thrsehold)
      Bivariate test (passing a Bonferroni thrsehold)

9. The bivariate results obtained with LAVA are also used by the script "phenogram_plot_input.R" for each combinations (AD and PD with and without proxies) to :
    
   -> Step 1 : Find the loci in the bivariate results that respect these thrsehold and save them in the appropriate directory. 
      - p<0.05 trehsold (save in /LAVA/phenoGram_plot/passing_bivariate_0.05 directory) 
      - Bonferroni thrsehold (save in /LAVA/phenoGram_plot/bivariate_bonferroni_thrsehold directory)
        
   -> Step 2 : Combine the loci into 1 file for each combinations of PD (sex-combined, sex-combined meta-analyzed, male and female) and AD (sex-combined) with and without proxies. 
      - use bivariate results to combined loci for PD and AD passing the univariate test saved (/LAVA/phenoGram_plot/passing univariate_test)
      - use loci passing bivariate test using p<0.05 (from step 1; /passing_bivariate_0.05 directory) to combine loci for PD and AD saved (/LAVA/phenoGram_plot/phenogram_data_005_threshold) 
      - use loci passing bivariate test using p<Bonferroni (from step 1; /bivariate_bonferroni_thrsehold directory) to combine loci for PD and AD saved (/LAVA/phenoGram_plot/bivariate_bonferroni_thrsehold)
      
11. The script "FUMA.R" combines the non-proxies PD/AD of all 6 kideny traits for all sexes and changes the SNPs to a gene set to be used in FUMA.
