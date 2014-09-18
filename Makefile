all:	GPL13534.bed

# Created on Mac terminal. Adapt for Linux, if needed

# Extract genomic coordinates of CpG sites on Illumina 450K HumanMethylation array
# Prerequisites: Use "Download full table..." button to get full annotation table from
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534. Gzip it
GPL13534.bed:	GPL13534-11288.txt.gz
				zcat < GPL13534-11288.txt.gz | sed '1,37d' | grep -e ^cg | cut -f1,12,13,17 | awk 'BEGIN {OFS="\t"} {if ($4 == "F") {S=$3; E=$3+50;} else {S=$3-50; E=$3;}} {print "chr"$2, S, E, $1}' | sort -k4 > GPL13534.bed