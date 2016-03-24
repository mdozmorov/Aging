all:	Illumina450K.bed Illumina27K.hg18.bed

# Created on Mac terminal. Adapt for Linux, if needed

Illumina450K.bed:	GPL13534-11288.txt.gz
	zcat < $< | sed '1,38d' | grep "^ch\|^cg" | cut -f1,12,13,17 | awk '{OFS="\t"} {if ($$4 == "F") {print "chr"$$2,$$3,$$3+49,$$1,".","+"} else {print "chr"$$2,$$3-48,$$3+1,$$1, ".","-"}}' | awk '{OFS="\t"} {if ($$2 <= 0) {print $$1,"0",$$3,$$4,$$5} else {print $0}}' | bedtools sort -i - | uniq > $@


Illumina27K.hg18.bed:	GPL8490-65.txt.gz
	zcat < $< | grep "^ch\|^cg" | sed '1,39d' | cut -f1,9,10,21 | awk '{OFS="\t"} {if ($$4 == "+") {print "chr"$$2,$$3,$$3+49,$$1,".","+"} else {print "chr"$$2,$$3-48,$$3+1,$$1, ".","-"}}' | awk '{OFS="\t"} {if ($$2 <= 0) {print $$1,"0",$$3,$$4,$$5} else {print $0}}' | bedtools sort -i - | uniq > $@
