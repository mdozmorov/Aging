all:	Horvath_cpg_all.bed

SHELL=/bin/bash

Horvath_cpg_pos.bed:	13059_2013_3156_MOESM3_ESM.csv
	join -t "	" -1 1 -2 4 <(sed '1,4d' $< | cut -d "," -f1-2 | awk '{FS=","} {OFS="\t"} $$2 > 0 {print $$1}' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Horvath_cpg_neg.bed:	13059_2013_3156_MOESM3_ESM.csv
	join -t "	" -1 1 -2 4 <(sed '1,4d' $< | cut -d "," -f1-2 | awk '{FS=","} {OFS="\t"} $$2 < 0 {print $$1}' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Horvath_cpg_all.bed:	Horvath_cpg_pos.bed Horvath_cpg_neg.bed
	cat $? > $@