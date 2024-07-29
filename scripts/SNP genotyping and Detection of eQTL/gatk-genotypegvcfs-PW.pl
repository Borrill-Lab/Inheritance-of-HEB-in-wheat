#!/usr/bin/perl -w

# Marek Glombik
# Adapted from Philippa Borrill

#
# Aim of the script is to run GATK haplotype caller on F5 RNAseq data

#### paths and references:
my $ref = "/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta";
my $input_dir = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/gatkout";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/gatkout";

### list of chrs to call variants on separately
my $list_of_chr = "/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/gatkout/intervals.list";

#############################

#open the list of chrs and go through the lines one by one to set off samtools and freebayes
open (INPUT_FILE, "$list_of_chr") || die "couldn't open the input file $list_of_chr!";
                    while (my $line = <INPUT_FILE>) {
                        chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";

#print "\nmy array: @array\n";
print "\narray element 1: @array[0]\n";


my $chr = $array[0];

chdir("$input_dir") or die "couldn't move to specific input directory $input_dir";


my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel freebayes tasks
#
#SBATCH -p nbi-medium
#SBATCH -t 2-00:00
#SBATCH --constraint=centos7
#SBATCH -c 4
#SBATCH --mem=100000
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/gatkout/gvcffinal-singlechr-gvcf-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWmap/picardout/gatkout/gvcffinal-singlechr-gvcf-%j.err
SLURM

 my $tmp_file = "$output_dir/gvcftmp/singlechrgatk.$chr";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $input_dir\n";

        print SLURM "set -e\n";
        print SLURM "source package 36919c64-b977-4ec7-b6d6-2bb35cb7808b\n";
        print SLURM "gatk GenomicsDBImport --genomicsdb-workspace-path finalPW-db_$chr -L $chr --batch-size 100 --sample-name-map gvcf_final_samplemap --tmp-dir $output_dir/gvcftmp -R $ref --reader-threads 4\n";
	print SLURM "gatk GenotypeGVCFs -R $ref -V gendb://finalPW-db_$chr -O all-finalPW-gvcf_$chr";

  close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;

}

            close(INPUT_FILE);
