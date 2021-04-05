#!/bin/bash

# Start in directory fastq
# Include info.tsv

#############################################################################################################################################

#Quality Check - Trimming

#for ((x = $1; x <= $2; x++))

#do
#	sample_fastq="0${x}.fastq.gz"
#	sample_trimmed="0${x}.T.fastq.gz"
#	fastq="0${x}.fastq.gz"

#	fastqc $sample_fastq

#	zcat $fastq | echo $((`wc -l`/4)) >> read_count.tmp

#	echo "Trimming - $sample_fastq.gz"

#	TrimmomaticSE -threads 5 $sample_fastq $sample_trimmed SLIDINGWINDOW:4:18 LEADING:28 TRAILING:28 MINLEN:36 >> log_Trimming 2>&1
#	paste info.tsv read_count.tmp > fastq_info.tsv

#done

#	wait

# Write number of reads in info_table
#paste info.tsv read_count.tmp > fastq_info.table.tsv
#rm read_count.tmp
#rm info.tsv

#mkdir trimmed
#mv *.T.fastq.gz trimmed/
#mv log_Trimming trimmed/

#mkdir base_quality
#mv *.html *.zip base_quality/

# SET DIR TRIMMED
#cd trimmed/

#fastqc *.gz

# SET DIR FASTQ
#cd ../

#############################################################################################################################################

# Mapping Hisat


#	echo "start mapping with hisat"
#	mkdir mapping
#	ln -s trimmed/* .

#for ((x = $1; x <= $2; x++))

#do
#	sample="0${x}.testing"
#	summary="summary_0${x}.txt"
#	fastq_input="0${x}.T.fastq.gz"
#	input_sam="0${x}.sam"
#	output_bam_tmp="0${x}.bam.tmp"
#	output_bam="0${x}.bam"
#	output_sorted_bam="0${x}.sorted.bam"

#echo "processing sample $fastq_input"

#hisat2 --threads 7 --summary-file $summary -x /media/dimbo/disk/GenOptics/data/genomes/mus_genome_hisat_index/mus_genome -U $fastq_input -S $input_sam

# keep only unique
#samtools view -h -F 4 $input_sam | awk 'substr($1, 0, 1)=="@" || $0 !~ /ZS:/' | samtools view -h -b > $output_bam_tmp

# remove dublicates
#samtools rmdup -s $output_bam_tmp $output_bam

# sort bam
#samtools sort $output_bam -o $output_sorted_bam

# index bam
#samtools index $output_sorted_bam

# remove tmp files
#rm $input_sam $output_bam.tmp $output_bam

#done

#	wait

#mv *.bam *.bai *.txt mapping/

#echo "mapping complete"

#############################################################################################################################################

# Create bar plots for uniquely and repeated mapped reads using R

# SET DIR MAPPING
#cd mapping/

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

# write R script for barplots and save it to file

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

	# BACK TO FASTQ
        cd ../

# Merge bam files (replicate 1 & 2) Treatment and control
# samtools merge Set8KO_TCP_input.merged.bam Set8KO_TCP_A_input.bam Set8KO_TCP_B_input.bam

#############################################################################################################################################

# QC Metrics (Deeptools) + Normalization

#mkdir qc_metrics
#cd qc_metrics/
#ln -s ../mapping/*.bam .
#ln -s ../mapping/*.bam.bai .

# CREATE MULTIBAM FILE
#multiBamSummary bins --bamfiles *.bam -p 7 -o multiBam.npz

# PLOT HEATMAP
#plotCorrelation -in multiBam.npz --corMethod spearman --colorMap Blues --plotTitle "Spearman Correlation of Read Counts" --whatToPlot heatmap --plotNumbers -o heatmap_SpearmanCorr_readCounts.png

# PLOT READ COVERAGE
#plotCoverage --bamfiles *.bam --skipZeros -p 7 --verbose -o Coverage_plot.png

# PLOT FINGERPRINT
#plotFingerprint -b *.bam -p 7 -plot finger_print_all.png

# CONVERT BAM TO BIGWIG (NORMALIZATION)
#mkdir normalization

#for ((x = $1; x <= $2; x++))

#do

#	sample_bam="0${x}.sorted.bam"
#	sample_bw="0${x}.bw"

#	echo "processing $sample_bam"

# Nomralize using BPM
#bamCoverage --bam $sample_bam -o $sample_bw --binSize 10 --outFileFormat bigwig --normalizeUsing BPM -p 7 --extendReads 200

#done

#wait

#mv *.bw normalization/

# Normalize using RPKM
# bamCoverage --bam $sample_bam -o $sample_bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPKM -p 7 --extendReads 200
# Normalize using RPGC
# bamCoverage --bam $sample_bam -o $sample_bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPGC -p 7 --effectiveGenomeSize 2407883318 --extendReads 200

# COMPUTEMATRIX

#cd normalization/

#	echo "Give input bigwick files"
#	echo "-----------------------"
#
#	echo "set treatment rep1"
#	read bw1
#	echo "set treatment rep2"
#	read bw2
#	echo "set control rep1"
#	read bw3
#	echo "set control rep2"
#	read bw4

#	output_matrix="matrix.gz"
#	output_regions="regions.bed"

#computeMatrix reference-point -b 5000 -a 5000 -R /media/dimbo/disk/GenOptics/data/genomes/genes_info/known_genes.mus.GENCODE.VM23.bed -S $bw1 $bw2 $bw3 $bw4 --skipZeros -o $output_matrix --outFileSortedRegions $output_regions -p 6

#	echo "compute matrix done"

# PLOT HEATMAP
#plotHeatmap -m $output_matrix  -out Heatmap.png

#mv matrix.gz regions.bed Heatmap.png ../

# RETURN TO QC_METRICS
cd ../

# RETURN TO FASTQ
cd ../

#############################################################################################################################################

# PEAK CALLING

#mkdir peak_calling

# call peak for TREATMENT (Set8KO) for replicate A and B
#macs2 callpeak --treatment Set8KO_TCP_A_H3K27ac.bam --control Set8KO_TCP_A_input.bam --nomodel --broad --format BAM --gsize mm -n rep1_H3K27ac_Set8Ko_A
#macs2 callpeak --treatment Set8KO_TCP_B_H3K27ac.bam --control Set8KO_TCP_B_input.bam --nomodel --broad --format BAM --gsize mm -n rep2_H3K27ac_Set8Ko_B

# call peak for INPUT (WTuntr) for replicate A and B
#macs2 callpeak --treatment WTuntr_A_H3K27ac.bam --control WTuntr_A_input.bam --nomodel --broad --format BAM --gsize mm -n rep1_H3K27ac_WTuntr_A
#macs2 callpeak --treatment WTuntr_B_H3K27ac.bam --control WTuntr_B_input.bam --nomodel --broad --format BAM --gsize mm -n rep2_H3K27ac_WTuntr_B

# Merge Peaks from Replicates
#bedtools intersect -a rep1_H3K27ac_Set8Ko_A_peaks.broadPeak -b rep2_H3K27ac_Set8Ko_B_peaks.broadPeak -wa > rep_H3K27ac_Set8KO_overlap_peaks.bed (H3K27ac_Set8KO_mergedPeaks.bed)
#bedtools intersect -a rep1_H3K27ac_WTuntr_A_peaks.broadPeak -b rep2_H3K27ac_WTuntr_B_peaks.broadPeak -wa > rep_H3K27ac_WTuntr_overlap_peaks.bed

# Exculde black listed regions
#bedtools intersect -v -a rep_H3K27ac_Set8KO_overlap_peaks.bed -b mm10-blacklist.v2.bed > rep_Filtered_H3K27ac_Set8KO_overlap_peaks.bed (H3K27ac_Set8KO_F.mergedPeaks.bed)
#bedtools intersect -v -a rep_H3K27ac_WTuntr_overlap_peaks.bed -b mm10-blacklist.v2.bed > rep_Filtered_H3K27ac_WTuntr_overlap_peaks.bed

# sort by p_values
#sort -k8 -g H3K27ac_Set8KO_FU.mergedPeaks.bed > H3K27ac_Set8KO_FU.mergedPeaks.sorted_P-values.bed

# Extract all p_values. Convert P-value (float numbers to integer)
#sort -k8 -g rep1_H3K27ac_Set8Ko_A_peaks.broadPeak | awk '{printf("%.0f\n",$8)}' > p_value.tmp

# replace new p_values
#awk 'BEGIN {OFS="\t"}; FNR==NR{a[NR]=$1;next}{$8=a[FNR]}1' p_value.tmp H3K27ac_Set8KO_FU.mergedPeaks.sorted_P-values.bed > H3K27ac_Set8KO_FU.mergedPeaks.sorted_P-values_FINAL.bed

# Remove unique IDS (in order to run ROSE)
#sort -u -k4 rep_Filtered_H3K27ac_WTuntr_overlap_peaks.bed > rep_Filtered_Unique_H3K27ac_WTuntr_overlap_peaks.bed (H3K27ac_Set8KO_FU.mergedPeaks.bed)

# convert to ROSE gff Format
#awk '{print $1,$4,$10,$2,$3,$10,$6,$10,$4}' OFS='\t' rep_H3K27ac_Set8KO.overlaps.broadPeak > rep_H3K27ac_Set8KO.overlaps.ROSE.gff (H3K27ac_Set8KO_FU.mergedPeaks.gff)

# Run ROSE. Identify SUPER ENHANCERS
#python ROSE_main.py --genome MM10 --i H3K27ac_Set8KO_FU.mergedPeaks.gff --rankby Set8KO_TCP_H3K27ac.merged.bam --control Set8KO_TCP_input.merged.bam --out SE_Set8KO/
#python ROSE_main.py --genome MM10 --i H3K27ac_WTuntr_FU.mergedPeaks.gff --rankby WTuntr_H3K27ac.merged.bam --control WTuntr_input.merged.bam --out SE_WTuntr/




