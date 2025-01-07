# Inheritance-of-HEB-in-wheat

<img src="https://github.com/user-attachments/assets/2e96c745-c24f-4d6f-bcf8-b3638a70438a" width="448" align="right">



## Background:

This repository contains scripts and links to data tables and raw sequencing data to a pre-print: "Rapid reprogramming and stabilisation of homoeolog expression bias in hexaploid wheat biparental populations" [(link)](https://www.biorxiv.org/content/10.1101/2024.08.01.606180v2.abstract)

In this study, we analyse Homoeolog Expression Bias inheritance patterns in two F<sub>5</sub> biparental populations with a common parent: Paragon (W10074 PF-1) x Charger (WPxCHA – 1001), and Paragon (W10074 PF-1) x Watkins (WATDE0228).


## Data availability:

All RNA-seq data were uploaded to [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/bioproject/?term=(PRJNA1128551)%20AND%20bioproject_sra[filter]%20NOT%20bioproject_gap[filter]) under the BioProject ID PRJNA1128551.

The wheat genome reference used to map/pseudoalign the RNA-seq data to is The IWGSC RefSeq v1.0 assembly available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/). 

The wheat transcriptome reference and genome annotation used to map SNPs to is The IWGSC RefSeq v1.1 annotation available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/).

The triad identification file "homoeologs_1_1_1_synt_and_non_synt.csv" is adapted from [HCtriads.tsv](https://opendata.earlham.ac.uk/wheat/under_license/toronto/Ramirez-Gonzalez_etal_2018-06025-Transcriptome-Landscape/data/TablesForExploration/) published as a part of dataset in [Ramírez-González et al. (2018)](https://www.science.org/doi/full/10.1126/science.aar6089).

All data tables necessary to run scripts associated with this project are located in the [figshare](https://figshare.com/projects/Reprogramming_and_stabilisation_of_homoeolog_expression_bias_in_hexaploid_wheat_biparental_populations/214495) repository, organised by Methods sections just like scripts (see below).

## Methodology:

Scripts are organised by Methods sections as outlined in the manuscript in the [scripts](scripts/) folder.


  ### 1. [HEB analysis](scripts/HEB/)


  ### 2. [Gene ontology term enrichment](scripts/Gene%20ontology%20term%20enrichment/)


  ### 3. [Variant calling, genotyping and association between genotype and homoeolog expression](scripts/Variant%20calling,%20genotyping%20and%20association%20between%20genotype%20and%20homoeolog%20expression/)

  
  ### 4. [SNP genotyping of F<sub>5</sub> lines for eQTL mapping](scripts/SNP%20genotyping%20and%20Detection%20of%20eQTL/)


  ### 5. [Quantification of gene expression levels and correction for batch effects](scripts/Quantification%20of%20gene%20expression%20levels%20and%20correction%20for%20batch%20effects/)


  ### 6. [Detection of eQTL](scripts/SNP%20genotyping%20and%20Detection%20of%20eQTL/)

_Note: Sections 4 and 6 are coupled together due to overlapping code sections._

## Figures:

The code for recreating figure panels is in the [figures](figures/) folder. The ready to plot data are located in the [figshare](https://figshare.com/projects/Reprogramming_and_stabilisation_of_homoeolog_expression_bias_in_hexaploid_wheat_biparental_populations/214495) repository associated with this project.
