#!/bin/bash

# Start in directory fastq
# Include info.tsv

#############################################################################################################################################

#	QUALITY BASE CHECK - TRIMMING

for ((x = $1; x <= $2; x++))

do
	sample_fastq="${x}.fastq.gz"
	sample_trimmed="${x}.T.fastq.gz"
	fastq="${x}.fastq.gz"

	fastqc $sample_fastq

	zcat $fastq | echo $((`wc -l`/4)) >> read_count.tmp

	echo "Trimming - $sample_fastq.gz"

	TrimmomaticSE -threads 5 $sample_fastq $sample_trimmed SLIDINGWINDOW:4:18 LEADING:28 TRAILING:28 MINLEN:36 >> log_Trimming 2>&1
	paste info.tsv read_count.tmp > fastq_info.tsv

done

	wait

#	Write number of reads in info_table
paste info.tsv read_count.tmp > fastq_info.table.tsv
rm read_count.tmp
rm info.tsv

mkdir trimmed
mv *.T.fastq.gz trimmed/
mv log_Trimming trimmed/

mkdir base_quality
mv *.html *.zip base_quality/

#	SET DIR TRIMMED
	cd trimmed/

fastqc *.gz

#	SET DIR FASTQ
	cd ../

#############################################################################################################################################

#	MAPPING HISAT


	echo "start mapping with hisat"
	echo "set path to reference genome"

	read genome

	mkdir mapping
	ln -s trimmed/* .

for ((x = $1; x <= $2; x++))

do
	sample="${x}.testing"
	summary="summary_${x}.txt"
	fastq_input="${x}.T.fastq.gz"
	input_sam="${x}.sam"
	output_bam_tmp="${x}.bam.tmp"
	output_bam="${x}.bam"
	output_sorted_bam="${x}.sorted.bam"

echo "processing sample $fastq_input"

# mapping Hisat2
hisat2 --threads 7 --summary-file $summary -x $genome -U $fastq_input -S $input_sam

# keep only unique
samtools view -h -F 4 $input_sam | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > $output_bam_tmp

# remove dublicates
samtools rmdup -s $output_bam_tmp $output_bam

# sort bam
samtools sort $output_bam -o $output_sorted_bam

# index bam
samtools index $output_sorted_bam

# remove tmp files
rm $input_sam $output_bam.tmp $output_bam

done

	wait

mv *.bam *.bai *.txt mapping/

echo "mapping complete"

#############################################################################################################################################

#	CREATE BAR PLOTS (R) FOR UNIQUE-REPEAT ALIGNMENTS

#	SET DIR MAPPING
	cd mapping/

for ((x = $1; x <= $2; x++))

do
	sum_file="summary_${x}.txt"

	uniq="uniq.tmp"
	repeat="repeat.tmp"
	uniq_repeat="uniq_repeat.tsv"
	identifier="identifier.tsv"
	id="${x}.sorted.bam"

	awk 'NR % 4 == 0' $sum_file | awk '{print $2}' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' >> $uniq
	awk 'NR % 5 == 0' $sum_file | awk '{print $2}' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?' >> $repeat
	echo $id | cut -f1 -d "." >> $identifier

done

	wait

	cat $uniq $repeat > $uniq_repeat
	ln -s ../fastq_info.table.tsv .
	rm *tmp

	num_samples=$(wc -l identifier.tsv | awk '{print $1}')

	awk '{print $14,$15}' fastq_info.table.tsv > identifier_new.tmp
	head -$num_samples identifier_new.tmp | tr -d "[:blank:]" > final_IDS.tsv

	working_directory=$(pwd)
	ncol=$(cat $uniq_repeat | echo $((`wc -l`/2)))

	paste identifier.tsv final_IDS.tsv > sample_IDs.tsv

#	Write R script for barplots and save it to file

	echo	"setwd(\"$working_directory\") 
		 colors = c(\"paleturquoise3\",\"red4\")
		 names = scan(\"identifier.tsv\", character(), quote = \"\")
		 values = scan(\"uniq_repeat.tsv\")
		 values = matrix(values, nrow = 2, ncol = $ncol, byrow = TRUE)
		 pdf(file = \"$working_directory/plot.pdf\", width = 7, height = 7)
		 barplot(values, names.arg = names, ylab = \"Alignment % Rate \", col = colors, density = 300, ylim = c(0,100) ,cex.lab = 1.4, las=3, cex.names = 1)
		 dev.off()" > R_Plot.R

		chmod 755 R_Plot.R
		Rscript R_Plot.R

	rm final_IDS.tsv identifier* uniq_repeat.tsv

#	SET DIR FASTQ
cd ../

# Merge bam files (replicate 1 & 2) Treatment and control
# samtools merge Set8KO_TCP_input.merged.bam Set8KO_TCP_A_input.bam Set8KO_TCP_B_input.bam

	echo "---------------------------------------------------------------------------------------------------------------"
	echo "----------------BAR PLOT FOR MAPPING RESULST COMPLETE!plot mapped results finished!----------------------------"
	echo "---------------------------------------------------------------------------------------------------------------"

#############################################################################################################################################

# QC Metrics (Deeptools) + Normalization

#############################################################################################################################################

# PEAK CALLING

mkdir peak_calling

#	SET DIR PEAK_CALLING
cd peak_calling/
ln -s /media/dimbo/disk/GenOptics/data/genomes/genes_info/mm10-blacklist.v2.bed .
ln -s ../mapping/*.bam .
ln -s ../mapping/*.bam.bai .

	echo "--------------------------------------------------------------------------------------------------------"
	echo "                                      CALL PEAKS USING MACS2                                            "
	echo "--------------------------------------------------------------------------------------------------------"

		while [[ $treatment != "none" ]]

		do

		echo "TYPE TREATMENT BAM FILE. else type [none]"
		read treatment

		echo "TYPE INPUT BAM FILE. else type [none]"
		read input

		echo "TYPE OUT PEAKS FILE. FORMAT:``condition_rep_A.`` else type [none]"
		read out

			if	[[ $treatment = "none" ]]
				then
				break
			else

			macs2 callpeak --treatment $treatment --control $input --nomodel --broad --broad-cutoff=0.01 --format BAM --gsize mm -n $out

		fi

		done


#		MERGE PEAKS FROM REPLICATES USING BEDTOOLS INTERSECT

		mkdir merged_peaks

		while [[ $repA != "none" ]]

		do

		echo "----------------------------------"
		echo "MERGE PEAKS FROM REPLICATE A AND B"
		echo "----------------------------------"

		echo "type peaks from replicate A. else type [none]"
		read repA

		echo "type peaks from replicate B. else type [none]"
		read repB

		echo "type overlapped peaks out file.  Format:  [overlapped_condition.tmp]"
		read out_file

			if      [[ $repA = "none" ]]
				then
				break
		else

		bedtools intersect -a $repA -b $repB -wa > $out_file



	# Exculde black listed regions

	final_peaks="$out_file.bed"

	echo "set path to blacklist_regions"
	read blacklist

	bedtools intersect -v -a $out_file -b $blacklist > $final_peaks

		fi

		done

	mv *.bed merged_peaks/
rm	*.tmp


#############################################################################################################################################

# IDENTIFY SUPER ENHANCERS (ROSE)

#############################################################################################################################################



#	SET DIR FASTQ
#	cd ../
#	CREATE SE DIR
#	mkdir SE
#	SET SE DIR
#	cd SE/

#	ln -s ../peak_calling/*.bed .

#	Remove unique IDS (in order to run ROSE)
#	sort -u -k4 rep_Filtered_H3K27ac_WTuntr_overlap_peaks.bed > rep_Filtered_Unique_H3K27ac_WTuntr_overlap_peaks.bed (H3K27ac_Set8KO_FU.mergedPeaks.bed)

# convert to ROSE gff Format
#awk '{print $1,$4,$10,$2,$3,$10,$6,$10,$4}' OFS='\t' rep_H3K27ac_Set8KO.overlaps.broadPeak > rep_H3K27ac_Set8KO.overlaps.ROSE.gff (H3K27ac_Set8KO_FU.mergedPeaks.gff)

# Run ROSE. Identify SUPER ENHANCERS
#python ROSE_main.py --genome MM10 --i H3K27ac_Set8KO_FU.mergedPeaks.gff --rankby Set8KO_TCP_H3K27ac.merged.bam --control Set8KO_TCP_input.merged.bam --out SE_Set8KO/
#python ROSE_main.py --genome MM10 --i H3K27ac_WTuntr_FU.mergedPeaks.gff --rankby WTuntr_H3K27ac.merged.bam --control WTuntr_input.merged.bam --out SE_WTuntr/




