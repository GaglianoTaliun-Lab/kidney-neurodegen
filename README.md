# kidney-neurodegen

# Global genetic correlation analysis (LDSC)

1. GWAS summary statistics files found in the "/sumstats_data" directory are formatted in a format acceptable by LDSC (ie: have columns "CHR" and "BP to be able to annotate rsids) using scripts "formatting_GWASsumstats.R" and "VCF_to_sumstats_format.R" and saves them in the "/processed_sumstats" directory. 
2. "Global_Genetic_Correlation_LDSC.sh" goes through all the different combinations using the new formatted GWAS summary statistics (/processed_sumstats directory) 


1. formatting_GWASsumstats.R and VCF_to_sumstats_format.R
   -> Takes original GWAS summary statistics files to change in the right format acceptable by LDSC : Have columns "CHR" and "BP" to be able to annotate rsids.
2. Global_Genetic_Correlation_LDSC.sh
   -> 
4. run_heritability.sh
   -> 
   
# Local genetic correlation analysis (LAVA)

1. input_info.R
   -> Uses a file named "case_control_input_meta_analyzed.csv" that has phenotype name, cases and controls of all AD/PD, PD meta (with and without proxies) in order to find the         paths and give "input_info_PD_meta.txt"
2. sample_overlap.R
   -> Step 1 : Uses all the results obtained from LDSC (in ldscore_regression directory) and combine them into 1 file and gives "ldsc_correlations_for_sample_overlap.txt" saved in       the directory "ldsc_corr".
   -> Step 2  : Creates the sample_overlap matrix for each correlation using "ldsc_correlations_for_sample_overlap.txt" (step 1 output) and saves them into the 
      "sample_overlap_files" directory.
3. 
