#!/usr/bin/perl -w

# Marek Glombik
# Adapted from Philippa Borrill

#
# Aim of script is to run picard and gatk for multiple samples for f5lines data

#### paths and references:
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout";
my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/";
my $input_for_picard = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picard_input_file.txt";


#######OUTPUT DIRECTORY, TMP DIRECTORY, SLURM_OUTPUT DIRECTORY MUST BE CREATED BEFORE RUNNING THE SCRIPT


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_picard") || die "couldn't open the input file $input_for_picard!";
        while (my $line = <INPUT_FILE>) {
      chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";
  
#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";
  
my $sample = $array[0];
  
  
chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";
  
  
my $SLURM_header = <<"SLURM";
#!/bin/bash
#SBATCH -p nbi-medium
#SBATCH -t 2-00:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --constraint=centos7
#SBATCH -J picard-gatk
#SBATCH --mail-type=none
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/picardslurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/picardslurm-%j.err
SLURM

 my $tmp_file = "$output_dir/tmppicard/picardgatk.$sample";
 
  
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;  
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";
  
  print SLURM "set -e\n";


  print SLURM "source package 7775013e-6118-4b30-86dd-0d7056a466fd\n";
  print SLURM "source package 36919c64-b977-4ec7-b6d6-2bb35cb7808b\n";


  print SLURM "picard MarkDuplicates -I $sample".".sorted.markdup.bam -O $output_dir/$sample"."_marked.bam -M $output_dir/$sample"."_metrics.txt\n";
  print SLURM "picard BuildBamIndex -I $output_dir/$sample"."_marked.bam\n";
  print SLURM "gatk SplitNCigarReads -R /jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta -I $output_dir/$sample"."_marked.bam -O  $output_dir/$sample"."_split.bam\n";
  print SLURM "picard AddOrReplaceReadGroups -I $output_dir/$sample"."_split.bam -O $output_dir/$sample"."_readgroup.bam -LB species -PL illumina -PU 1 -SM $output_dir/$sample\n";
  print SLURM "picard BuildBamIndex -I $output_dir/$sample"."_readgroup.bam\n";

    close SLURM;
    system("sbatch $tmp_file");
 # unlink $tmp_file;
  
}

close(INPUT_FILE);
