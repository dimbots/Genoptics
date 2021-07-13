#!/bin/bash
# QC Metrics (Deeptools) + Normalization

#	START IN FASTQ DIRECTORY
#	SET DIR QC_METRICS
	mkdir qc_metrics
	cd qc_metrics/

	ln -s ../mapping/*.bam .
	ln -s ../mapping/*.bam.bai .

	tput setaf 1; tput bold; echo "RENAME FILES FROM 0869... TO CONDITIONS (e.g Set8KO_A_input.) IN DIRECTORY (qc_metrics). WHEN COMPLETE TYPE [done]."
	read response

# CREATE MULTIBAM FILE
multiBamSummary bins --bamfiles *.bam -p 7 -o multiBam.npz

# PLOT HEATMAP
plotCorrelation -in multiBam.npz --corMethod spearman --colorMap Blues --plotHeight 11.5 --plotWidth 13  --whatToPlot heatmap --plotNumbers -o SpearmanCor_readCounts.png

# PLOT READ COVERAGE
plotCoverage --bamfiles *.bam --skipZeros -p 7 --verbose -o Coverage_plot.png

# PLOT FINGERPRINT
plotFingerprint -b *.bam -p 7 -plot finger_print_all.png

#########################################################################################

#	NORMALIZE BAM FILES

#	CONVERT BAM TO BIGWIG (NORMALIZATION)



	tput setaf 1; tput bold; echo "--------------------------------------------------------------------------"
	tput setaf 1; tput bold; echo "				START NORMALIZATION				"
	tput setaf 1; tput bold; echo "--------------------------------------------------------------------------"


	tput setaf 2; tput bold; echo "TYPE METHOD OF NORMALIZATION. RPKM / BPM / RPGC"
	read method

	if
	[[ $method = "BPM" ]]
			then

			for i in $(ls *.bam)

			do
				bw="$i.bw"
				bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing BPM -p 7 --extendReads 200
			done
	fi



		if

		[[ $method = "RPKM" ]]
		then
			for i in $(ls *.bam)

			do

			bw="$i.bw"
			bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPKM -p 7 --extendReads 200

			done
		fi



	if
	[[ $method = "RPGC" ]]
	then

		for i in $(ls *.bam)

		do
		bw="$i.bw"
		bamCoverage --bam $i -o $bw --binSize 10 --outFileFormat bigwig --normalizeUsing RPGC -p 7 --effectiveGenomeSize 2407883318 --extendReads 200

		done
	fi

mkdir normalization
mv *.bw normalization/

#########################################################################################

# 			COMPUTEMATRIX

#	SET DIR NORMALIZATION
	cd normalization/

	tput setaf 1; tput bold; echo "          START COMPUTEMATRIX            "

	tput setaf 1; tput bold; echo "-----------------------------------------"
	tput setaf 1; tput bold; echo "   TYPE BIGWIG FILES (TREATMENT-INPUT)"
	tput setaf 1; tput bold; echo "-----------------------------------------"

	tput setaf 2; tput bold; echo "TYPE TREATMENT BW FILE. ELSE TYPE [none]"
	read bw1
	tput setaf 2; tput bold; echo "TYPE INPUT BW FILE. ELSE TYPE [none]"
	read bw2
	tput setaf 2; tput bold; echo "TYPE TREATMENT BW FILE. ELSE TYPE [none]"
	read bw3
	tput setaf 2; tput bold; echo "TYPE INPUT BW FILE. ELSE TYPE [none]"
	read bw4
	tput setaf 2; tput bold; echo "TYPE OUT MATRIX FILE. E.G (matrix.0865_0866)"
	read out_matrix
	tput setaf 2; tput bold; echo "TYPE OUT REGIONG FILE. E.G (out_regions.0865_0866)"
	read out_regions

	tput setaf 2; tput bold; echo "Set path to genes.bed file"
	read genes

	tput setaf 2; tput bold; echo "TYPE NUMBER OF BASE PAIRS BEFORE TSS"
	read before_tss
	tput setaf 2; tput bold; echo "TYPE NUMBER OF BASE PAIRS AFTER TSS"
	read after_tss

	if [[ "$bw3" != "none" && "bw4" != "none" ]]
		then
	tput setaf 1; tput bold; echo "RUNNING WITH 4 SAMPLES"
		computeMatrix reference-point -b $before_tss -a $after_tss -R $genes -S $bw1 $bw2 $bw3 $bw4 --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 7

		out_plot="Heatmap_$out_matrix"
		plotHeatmap -m $out_matrix -out $out_plot --colorMap Blues
	else
	tput setaf 1; tput bold; echo "RUNNING WITH 2 SAMPLES"
		computeMatrix reference-point -b $before_tss -a $after_tss -R $genes -S $bw1 $bw2 --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 7

		out_plot="Heatmap_$out_matrix"
		plotHeatmap -m $out_matrix -out $out_plot --colorMap Blues
	fi

mv $out_plot ../







