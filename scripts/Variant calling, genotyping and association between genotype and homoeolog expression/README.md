# Variant calling, genotyping and association between genotype and homoeolog expression

Data tables necessary to run the scripts are available on figshare

## Scripts:

trim.pl - performs trimming of the RNA-seq data

hisat.pl - performs mapping of trimmed RNa-seq data to the genome reference

picard-markdup-index.pl - performs preparation of bamfiles for variant calling

haplotype-caller.pl - performs haplotype calling on individiual samples

gatk-genotypegvcfs.pl - performs joint genotyping on samples

snp_genotyping_triads.R - associates genotype information with triads/homoeologs

association_genotype_expression.R - analyses the association between homoeolog expression and inherited genotype
