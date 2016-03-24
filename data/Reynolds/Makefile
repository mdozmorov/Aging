all:	Reynolds_cpg_tcells_all.bed Reynolds_cpg_monocytes_all.bed

SHELL=/bin/bash

Reynolds_cpg_tcells_pos.bed:	ncomms6366-s2.txt
	join -t "	" -1 1 -2 4 <(cat $< | sed '1,4d' | cut -f1,6 | awk '{if ( $$2 > 0 ) { print $$1 } }' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Reynolds_cpg_tcells_neg.bed:	ncomms6366-s2.txt
	join -t "	" -1 1 -2 4 <(cat $< | sed '1,4d' | cut -f1,6 | awk '{if ( $$2 < 0 ) { print $$1 } }' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Reynolds_cpg_monocytes_pos.bed:	ncomms6366-s3.txt
	join -t "	" -1 1 -2 4 <(cat $< | sed '1,4d' | cut -f1,9 | awk '{if ( $$2 > 0 ) { print $$1 } }' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Reynolds_cpg_monocytes_neg.bed:	ncomms6366-s3.txt
	join -t "	" -1 1 -2 4 <(cat $< | sed '1,4d' | cut -f1,9 | awk '{if ( $$2 < 0 ) { print $$1 } }' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Reynolds_cpg_tcells_all.bed:	Reynolds_cpg_tcells_pos.bed Reynolds_cpg_tcells_neg.bed
	cat $? > $@

Reynolds_cpg_monocytes_all.bed:	Reynolds_cpg_monocytes_pos.bed Reynolds_cpg_monocytes_neg.bed
	cat $? > $@

