#!/usr/bin/perl -w

# Marek Glombik
# Adapted from Philippa Borrill

#
# Aim of script is to run hisat2 for multiple samples for f5lines data

#### paths and references:
my $path = '/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1';
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Arun/f6_lines/mapping/fastq_trim/";
my $ref = "$path/161010_Chinese_Spring_v1.0_pseudomolecules_parts";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Marek/mg_mapping/hisat2/";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $output_dir):
my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Arun/f6_lines/mapping/fastq_trim/";
my $input_for_hisat = "/jic/scratch/groups/Philippa-Borrill/Arun/f6_lines/mapping/fastq_trim/for_hisat_sample_list.txt";


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_hisat") || die "couldn't open the input file $input_for_hisat!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";
  
#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";
  
my $sample = $array[0];
my $f1 = $array[1];
my $f2 = $array[2];
  
  
chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";
  
  
my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel hisat2 tasks
#
#SBATCH -p jic-medium
#SBATCH -t 2-00:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH -J hisat2sep
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Marek/mg_mapping/hisat2/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Marek/mg_mapping/hisat2/slurm_output/%x.%N.%j.err
SLURM
  
 my $tmp_file = "$output_dir/tmp/hisat2.$sample";
  
  
  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";
  
  print SLURM "set -e\n";
  
  print SLURM "source package f9c1e0c5-d0e8-4ba0-9edd-88235400fa13\n";
  
    ### this part all needs editing!!!! ###
    print SLURM "hisat2 --rg-id $sample --rg 'SM:$sample' --rg 'LB:$sample' -p 4  --dta -x $ref -1 $f1 -2 $f2 -S $output_dir/$sample".".sam\n";
    print SLURM "source package aeee87c4-1923-4732-aca2-f2aff23580cc\n";
    print SLURM "samtools sort -@ 4 -n -T $output_dir/$sample"."temp_namesort.bam -O bam -o $output_dir/$sample"."sorted.bam $output_dir/$sample".".sam\n";
    print SLURM "rm $output_dir/$sample".".sam\n";
    print SLURM "samtools fixmate -m $output_dir/$sample"."sorted.bam $output_dir/$sample"."fixmate.bam\n";
    print SLURM "rm $output_dir/$sample"."sorted.bam\n";
        print SLURM "samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b -o $output_dir/$sample"."q20.bam  $output_dir/$sample"."fixmate.bam\n";
    print SLURM "rm $output_dir/$sample"."fixmate.bam\n";
    print SLURM "samtools sort -@ 4 -T $output_dir/$sample"."temp_pos_sort.bam -O bam -o $output_dir/$sample"."pos_sort.bam $output_dir/$sample"."q20.bam\n";
        print SLURM "rm $output_dir/$sample"."q20.bam\n";
    print SLURM "samtools markdup -r $output_dir/$sample"."pos_sort.bam -T $output_dir/$sample"."temp_markdup.bam $output_dir/$sample".".sorted.markdup.bam\n";
    print SLURM "rm $output_dir/$sample"."pos_sort.bam\n";
    print SLURM "samtools index $output_dir/$sample".".sorted.markdup.bam\n";
  
  close SLURM;
    system("sbatch $tmp_file");
  # unlink $tmp_file;
  
}

close(INPUT_FILE);
