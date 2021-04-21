#!/bin/bash

# Start in directory fastq
# Include info.tsv

#############################################################################################################################################

#	QUALITY BASE CHECK - TRIMMING

for ((x = $1; x <= $2; x++))

do
	sample_fastq="0${x}.fastq.gz"
	sample_trimmed="0${x}.T.fastq.gz"
	fastq="0${x}.fastq.gz"

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
	mkdir mapping
	ln -s trimmed/* .

for ((x = $1; x <= $2; x++))

do
	sample="0${x}.testing"
	summary="summary_0${x}.txt"
	fastq_input="0${x}.T.fastq.gz"
	input_sam="0${x}.sam"
	output_bam_tmp="0${x}.bam.tmp"
	output_bam="0${x}.bam"
	output_sorted_bam="0${x}.sorted.bam"

echo "processing sample $fastq_input"

hisat2 --threads 7 --summary-file $summary -x /media/dimbo/disk/GenOptics/data/genomes/mus_genome_hisat_index/mus_genome -U $fastq_input -S $input_sam

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
	sum_file="summary_0${x}.txt"

	uniq="uniq.tmp"
	repeat="repeat.tmp"
	uniq_repeat="uniq_repeat.tsv"
	identifier="identifier.tsv"
	id="0${x}.sorted.bam"

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
	bedtools intersect -v -a $out_file -b mm10-blacklist.v2.bed > $final_peaks

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




#!/bin/bash
# QC Metrics (Deeptools) + Normalization

#	SET DIR QC_METRICS
#	mkdir qc_metrics
	cd qc_metrics/

#ln -s ../mapping/*.bam .
#ln -s ../mapping/*.bam.bai .

#	echo "RENAME FILES FROM 0869... TO CONDITIONS (e.g Set8KO_A_input.) IN DIRECTORY (qc_metrics). WHEN COMPLETE TYPE [done]."
#	read response

# CREATE MULTIBAM FILE
#multiBamSummary bins --bamfiles *. -p 7 -o multiBam.npz

# PLOT HEATMAP
#plotCorrelation -in multiBam.npz --corMethod spearman --colorMap Blues --plotHeight 11.5 --plotWidth 13  --whatToPlot heatmap --plotNumbers -o SpearmanCor_readCounts.png

# PLOT READ COVERAGE
#plotCoverage --bamfiles *. --skipZeros -p 7 --verbose -o Coverage_plot.png

# PLOT FINGERPRINT
#plotFingerprint -b *. -p 7 -plot finger_print_all.png

#########################################################################################

#	NORMALIZE BAM FILES

#	CONVERT BAM TO BIGWIG (NORMALIZATION)



#	echo "--------------------------------------------------------------------------"
#	echo "				START NORMALIZATION				"
#	echo "--------------------------------------------------------------------------"


#		echo "TYPE METHOD OF NORMALIZATION. RPKM / BPM / RPGC"
#		read method

#		if
#			[[ $method = "BPM" ]]
#			then

#				for i in $(ls *.)

#                			do
#					bw="$i.bw"
#					bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing BPM -p 7 --extendReads 200
#				done
#				fi

#		if
#			[[ $method = "RPKM" ]]
#			then

#				for i in $(ls *.)

#                                	do
 #                               	bw="$ibw"
#					bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPKM -p 7 --extendReads 200
#				done
#				fi

#		if
#			[[ $method = "RPGC" ]]
#			then
#
#				for i in $(ls *.)

#					do
#					bw="$ibw"
#					bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPGC -p 7 --effectiveGenomeSize 2407883318 --extendReads 200
#				done
#				fi

#mkdir normalization
#mv *.bw normalization/

#########################################################################################

# 			COMPUTEMATRIX

#	SET DIR NORMALIZATION
	cd normalization/

#	echo "          START COMPUTEMATRIX            " 

#	echo "-----------------------------------------"
#	echo "   TYPE BIGWIG FILES (TREATMENT-INPUT)"
#	echo "-----------------------------------------"

#		echo "TYPE TREATMENT BW FILE. ELSE TYPE [none]"
#		read bw1

#		echo "TYPE INPUT BW FILE. ELSE TYPE [none]"
#		read bw2

#		echo "TYPE TREATMENT BW FILE. ELSE TYPE [none]"
#		read bw3

#		echo "TYPE INPUT BW FILE. ELSE TYPE [none]"
#		read bw4

#		echo "TYPE OUT MATRIX FILE. E.G (matrix.0865_0866)"
#		read out_matrix

#		echo "TYPE OUT REGIONG FILE. E.G (out_regions.0865_0866)"
#		read out_regions


#		if [[ "$bw3" != "none" && "bw4" != "none" ]]
#			then
#			echo "RUNNING WITH 4 SAMPLES"
#			computeMatrix reference-point -b 10000 -a 10000 -R /media/dimbo/disk/GenOptics/data/genomes/genes_info/known_genes.mus.GENCODE.VM23.bed -S $bw1 $bw2 $bw3 $bw4 --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 7

			echo "NAMES IN PLOT"
			echo "TYPE SAMPLE A TREATMENT"
			read a
			echo "TYPE SAMPLE A INPUT"
			read b
			echo "SAMPLE B TREATMENT"
			read c
			echo "SAMPLE B INPUT"
			read d
			
			# THIS IS OPTIONAL DELETED AFTER TEST
			echo "TYPE MATRIX OUT FILE NAME"
			read $output_matrix

			out_plot="Heatmap_$out_matrix"
			plotHeatmap -m $output_matrix -out $out_plot --colorMap Blues --samplesLabel $a $b $c $d


#		else
#			echo "RUNNING WITH 2 SAMPLES"
#
#			echo "NAMES IN PLOT"
#			echo "TYPE SAMPLE A TREATMENT"
#			read a
#			echo "TYPE SAMPLE A INPUT"
#			read b

			# THIS IS OPTIONAL DELETE AFTER TEST
#			echo "TYPE MATRIX OUT FILE NAME"
#			read $output_matrix

#			out_plot="Heatmap_$out_matrix"
#			computeMatrix reference-point -b 10000 -a 10000 -R /media/dimbo/disk/GenOptics/data/genomes/genes_info/known_genes.mus.GENCODE.VM23.bed -S $bw1 $bw2 --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 7

#			echo "NAMES IN PLOT"
#			echo "TYPE SAMPLE A TREATMENT"
#			read a
#			echo "TYPE SAMPLE A INPUT"
#			read b

#			echo "TYPE MATRIX OUT FILE NAME"
#			read $output_matrix

#			out_plot="Heatmap_$out_matrix"
#			plotHeatmap -m $output_matrix -out $out_plot --colorMap Blues --samplesLabel $a $b



#		fi



# PLOT HEATMAP

#	out_plot="Heatmap_$out_matrix"

#	echo "SAMPLE A TREATMENT"
#	read a
#	echo "SAMPLE A INPUT"
#	read b
#	echo "SAMPLE B TREATMENT "
#	read c
#	echo "sample b input"
#	read d

#	plotHeatmap -m $output_matrix -out $out_plot --colorMap Blues --samplesLabel $a $b $c $d

#	mv $out_plot ../







