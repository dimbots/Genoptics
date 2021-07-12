#!/bin/bash


for ((x = $1; x <= $2; x++))

do

	input="${x}.sorted.bam"
	header="${x}.tmp"

	samtools view -H $input > $header

	shuffled="${x}.sam.tmp"
	
	echo "type the number of least mapped reads in the data set"
	read number

	# -n (this is the sample that has the least mapped reads. (extracted from bam file))
	samtools view $input | shuf | head -n $number > $shuffled

	unsorted="${x}.downsampled.tmp"

	cat $header $shuffled > $unsorted

	sorted="${x}.downsampled.bam"

	samtools sort $unsorted -o $sorted

	samtools index $sorted

	bw="${x}.bw"
	bamCoverage -b $sorted -o $bw 

	rm $shuffled $unsorted
done

	rm *.bam *.bai *.tmp
