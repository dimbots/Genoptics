bedtools makewindows -g mm10.chrom.sizes -w 5000 | awk '$1 == "chr2" {print $0}' > chr2.bed
