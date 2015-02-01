all:	Illumina450K.bed Illumina27K.hg18.bed

# Created on Mac terminal. Adapt for Linux, if needed

Illumina450K.bed:	GPL13534-11288.txt.gz
					zcat < GPL13534-11288.txt.gz | sed '1,37d' | cut -f1,7,12,17,34,35,17 | grep -e ^cg | awk 'BEGIN {OFS="\t"} {if ($$4 == "F") {S="+";} else {S="-";}} {print "chr"$$3,$$5,$$6,$$1,$$2,S;}' | sort -k4,4 | uniq > $@

Illumina27K.hg18.bed:	GPL8490-65.txt.gz
	zcat < $< | sed '1,39d' | cut -f1,9,10,21 | awk 'BEGIN {OFS="\t"} {print "chr"$$2,$$3,$$3+123,$$1,"0",$$4}' | sort -k4,4 | uniq > $@