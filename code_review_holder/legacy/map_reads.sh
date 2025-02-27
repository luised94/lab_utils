#!/usr/bin/env bash

# Written by JBK on 19.06.14
# Example usage:
#$ ./map_reads.sh reads/r1.fq reads/r2.fq genome/s288c.fa tmp/ pileup.txt stats.txt

# Warning: genome must be prepped before this. 
# Do this by running:
# ./prep_genome <genome_file>

# Get arguments
read1_file=$1 # Read1 fastq file
read2_file=$2 # Read2 fastq file
genome_file=$3 # Genome fasta file
tmp_prefix=$4 # Prefix for all tmp files
out_pileup_file=$5 # Output pileup file
out_stats_file=$6

# Map reads
echo "Mapping reads..." >&2

# Create temporary file prefix
prefix="$(mktemp -p $tmp_prefix)"
echo "Using tmp file prefix: $prefix" >&2

# Map paired-end reads
bwa mem $genome_file $read1_file $read2_file > $prefix.sam

# Convert text sam file to binary bam file
samtools view -S -b $prefix.sam > $prefix.bam

# Filter reads based on various criteria
samtools view -b -f 2 -F 524 $prefix.bam > $prefix.filtered.bam

# Fuse paired end reads in bed format 
bedtools bamtobed -bedpe -i $prefix.filtered.bam > $prefix.bed

# Sort bed file; required before genomcov is used
bedtools sort -i $prefix.bed > $prefix.sorted.bed

# Finally, compute coverage in bedgraph format and send to standard output
bedtools genomecov -bg -i $prefix.sorted.bed -g $genome_file.genome \
    > $out_pileup_file

# Write stats file
samtools flagstat $prefix.bam > $out_stats_file
# This isn't necessary, because only the "properly paired" reads are used
# after filtering
# samtools flagstat $prefix.filtered.bam >> $out_stats_file

# Clear tmp files
echo "Removing files with prefix: $prefix" >&2
#rm $prefix*

echo "Done." >&2