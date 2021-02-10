#!/bin/bash

for ((x = $1; x <= $2; x++))

do
	sample="0${x}.tesing"
	summary="summary_0${x}.txt"
	fastq_input="0${x}.fastq"
	input_sam="0${x}.sam"
	output_bam_tmp="0${x}.bam.tmp"
	output_bam="0${x}.bam"
	output_sorted_bam="0${x}.sorted.bam"

echo "testing"
echo $sample

# mapping with hisat2
hisat2 -p 7 --summary-file $summary -x ../../../../Genomes_Fasta/mus_Genome/hista2_index/mus_genome -U $fastq_input -S $input_sam

# keep only unique
samtools view -h -F 4 $input_sam | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > $output_bam_tmp

# remove dublicates
samtools rmdup -s $output_bam.tmp $output_bam

# sort bam
samtools sort $output_bam -o $output_sorted_bam

# index bam
samtools index $output_sorted_bam

# remove tmp files
rm $input_sam $output_bam.tmp $output_bam


done
