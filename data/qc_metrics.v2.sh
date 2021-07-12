#!/bin/bash
# QC Metrics (Deeptools) + Normalization

#	SET DIR QC_METRICS
#	mkdir qc_metrics
#	cd qc_metrics/

#ln -s ../mapping/*.bam .
#ln -s ../mapping/*.bam.bai .

#	echo "RENAME FILES FROM 0869... TO CONDITIONS (e.g Set8KO_A_input.) IN DIRECTORY (qc_metrics). WHEN COMPLETE TYPE [done]."
#	read response

# CREATE MULTIBAM FILE
#multiBamSummary bins --bamfiles *.bam -p 7 -o multiBam.npz

# PLOT HEATMAP
#plotCorrelation -in multiBam.npz --corMethod spearman --colorMap Blues --plotHeight 11.5 --plotWidth 13  --whatToPlot heatmap --plotNumbers -o SpearmanCor_readCounts.png

# PLOT READ COVERAGE
#plotCoverage --bamfiles *.bam --skipZeros -p 7 --verbose -o Coverage_plot.png

# PLOT FINGERPRINT
#plotFingerprint -b *.bam -p 7 -plot finger_print_all.png

#########################################################################################

#	NORMALIZE BAM FILES

#	CONVERT BAM TO BIGWIG (NORMALIZATION)



	echo "--------------------------------------------------------------------------"
	echo "				START NORMALIZATION				"
	echo "--------------------------------------------------------------------------"


	echo "TYPE METHOD OF NORMALIZATION. RPKM / BPM / RPGC"
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
#	cd normalization/

#	echo "          START COMPUTEMATRIX            "

#	echo "-----------------------------------------"
#	echo "   TYPE BIGWIG FILES (TREATMENT-INPUT)"
#	echo "-----------------------------------------"

#	echo "TYPE TREATMENT BW FILE. ELSE TYPE [none]"
#	read bw1
#	echo "TYPE INPUT BW FILE. ELSE TYPE [none]"
#	read bw2
#	echo "TYPE TREATMENT BW FILE. ELSE TYPE [none]"
#	read bw3
#	echo "TYPE INPUT BW FILE. ELSE TYPE [none]"
#	read bw4
#	echo "TYPE OUT MATRIX FILE. E.G (matrix.0865_0866)"
#	read out_matrix
#	echo "TYPE OUT REGIONG FILE. E.G (out_regions.0865_0866)"
#	read out_regions

#	echo "Set path to genes.bed file"
#	read genes


#	if [[ "$bw3" != "none" && "bw4" != "none" ]]
#		then
#		echo "RUNNING WITH 4 SAMPLES"
#		computeMatrix reference-point -b 5000 -a 5000 -R $genes -S $bw1 $bw2 $bw3 $bw4 --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 7

#		out_plot="Heatmap_$out_matrix"
#		plotHeatmap -m $out_matrix -out $out_plot --colorMap Blues
#	else
#		echo "RUNNING WITH 2 SAMPLES"
#		computeMatrix reference-point -b 5000 -a 5000 -R $genes -S $bw1 $bw2 --skipZeros -o $out_matrix --outFileSortedRegions $out_regions -p 7

#		out_plot="Heatmap_$out_matrix"
#		plotHeatmap -m $out_matrix -out $out_plot --colorMap Blues
#	fi

#mv $out_plot ../







