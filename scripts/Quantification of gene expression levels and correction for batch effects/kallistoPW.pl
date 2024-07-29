#!/usr/bin/perl -w

# Aim of script is to run kallisto on RNA-seq for multiple samples to a common reference to calculate expression levels
# Adapted from Philippa Borrill

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06';
my $ref = "$path/IWGSC_v1.1_ALL_20170706_transcripts.fasta";
my $index = "$path/IWGSC_v1.1_ALL_20170706_transcripts.fasta_index";

#############################



# NB make index by kallisto index  -i IWGSC_v1.1_ALL_20170706_transcripts_index IWGSC_v1.1_ALL_20170706_transcripts.fasta

my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PW/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWkallisto/";

### lists of samples (text file containing directory/subdirectory with .fastq 
my $input_for_kallisto = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PW/PWinput.txt";


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$read_path_triticum") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_kallisto") || die "couldn't open the input file $input_for_kallisto!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";

#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";

my $sample = $array[0];
my $fastq1 = $array[1];
my $fastq2 = $array[2];


chdir("$read_path_triticum") or die "couldn't move to specific read directory";


my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel hisat2 tasks
#
#SBATCH -p nbi-medium
#SBATCH -t 0-05:00
#SBATCH -c 8
#SBATCH --mem=30000
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWkallisto/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWkallisto/%x.%N.%j.err


SLURM

 my $tmp_file = "$output_dir/tmp/kallisto.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";


  print SLURM "set -e\n";

	print SLURM "source package /nbi/software/testing/bin/kallisto-0.46.1\n";
  print SLURM "source package aeee87c4-1923-4732-aca2-f2aff23580cc\n";

	print SLURM "kallisto quant -i $index -o $output_dir/$sample -t 8 --pseudobam $fastq1 $fastq2 \n";
	print SLURM "samtools sort -T $output_dir/$sample/$sample".".bam -O bam -o $output_dir/$sample".".sorted.bam -@ 8 $output_dir/$sample/pseudoalignments.bam\n";
  print SLURM "samtools index $output_dir/$sample".".sorted.bam $output_dir/$sample".".sorted.bai\n";
  print SLURM "rm $output_dir/$sample/pseudoalignments.bam\n";

	close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;

}

	    close(INPUT_FILE);
