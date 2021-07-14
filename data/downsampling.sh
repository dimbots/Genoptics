#!/bin/bash

# DOWNSAMPLING BASED ON SAMPLE WITH THE LESS NUMBER OF READS

# CREATE DIR DOWNSAMPLING
mkdir downsampling
# SET DIR TO DOWNSAMPLING
cd downsampling/

ln -s ../mapping/*.bam .

#	tput setaf 2; tput bold; echo "TUPE THE NUMBER OF LEAST MAPPED READS IN THE DATA SET"
#	read number

	cat ../mapping/summary_08* > mapped_reads.tmp
	number=$(awk 'NR%3==1' ../mapping/mapped_reads.tmp | grep -v "^ " | cut -d ' ' -f 1 | sort -n | head -1)

	tput setaf 2; tput bold; echo "MAPPED SAMPLE WITH LESS NUMBER OF READS IS: $number"

	tput setaf 1; tput bold; echo "------------------"
	tput setaf 1; tput bold; echo "START DOWNSAMPLING"
	tput setaf 1; tput bold; echo "------------------"

for ((x = $1; x <= $2; x++))

do

	input="0${x}.sorted.bam"
	header="0${x}.tmp"

	samtools view -H $input > $header

	shuffled="0${x}.sam.tmp"

	# -n (this is the sample that has the least mapped reads. (extracted from bam file))
	samtools view -@ 7 $input | shuf | head -n $number > $shuffled

	unsorted="0${x}.downsampled.tmp"

	cat $header $shuffled > $unsorted

	sorted="0${x}.downsampled.bam"

	samtools sort -@ 7 $unsorted -o $sorted

	samtools index -@ 7 $sorted

	bw="0${x}_downsampled.bw"
	bamCoverage -p 7 -b $sorted -o $bw 

	rm $shuffled $unsorted

done

	rm *.bam *.bai *.tmp
