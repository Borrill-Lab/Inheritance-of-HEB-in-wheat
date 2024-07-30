# HEB

Data tables necessary to run the scripts are available on [figshare](https://figshare.com/articles/dataset/Data_tables_for_HEB/26405026)

## Scripts:

kallisto.pl - run pseudoalignment and quantification on RNA-seq data

tximport_summarise_counts_tpm_per_gene.R - extract raw counts and TPM values

assign_homoeolog_expression_bias_categories.R - calculates HEB categories based on [Ramírez-González et al. (2018)](https://www.science.org/doi/full/10.1126/science.aar6089)


F5_lines_heb.R - calculates HEB as CV and initial bias distance calculations in F5 lines

F5_lines_HEB_additional_cv_analysis.R - calculates CV relationship to bias distance, position on chromosome, creates panels for figure 1-2

bias_dist_stats.R - performs bias distance calculations and creates panels for figure 2
