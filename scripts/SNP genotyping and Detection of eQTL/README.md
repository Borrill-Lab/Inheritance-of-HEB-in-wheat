# SNP genotyping of F<sub>5</sub> lines

RNA-seq data necessary to run the scripts are available on SRA under the BioProject ID [PRJNA1128551](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1128551/).

## Scripts:

trimPC/PW60.pl - performs trimming of the RNA-seq data

hisatPC/PW.pl - performs mapping of trimmed RNa-seq data to the genome reference

picard-markdup-index-PC/PW.pl - performs preparation of bamfiles for variant calling

haplotype-caller-PC/PW.pl - performs haplotype calling on individiual samples

gatk-genotypegvcfs.PC/PW.pl - performs joint genotyping on samples


# Detection of eQTL

Data tables necessary to run the scripts are available on [figshare](https://figshare.com/account/projects/214495/articles/28308998?file=52023947)

## Scripts:

LMM_eQTL_PC_PW.R - runs the SNP filtering and subsequent LMM and eQTL analysis
