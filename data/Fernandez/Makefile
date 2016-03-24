all:	Fernandez_MSCs_cpg_all.bed Fernandez_blood_cpg_all.bed

SHELL=/bin/bash

Fernandez_MSCs_cpg_pos.bed:	Supplemental_Table_2.txt
	join -t "	" -1 1 -2 4 <(sed '1,2d' $< | cut -f1 | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Fernandez_MSCs_cpg_neg.bed:	Supplemental_Table_3.txt
	join -t "	" -1 1 -2 4 <(sed '1,2d' $< | cut -f1 | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Fernandez_MSCs_cpg_all.bed:	Fernandez_MSCs_cpg_pos.bed Fernandez_MSCs_cpg_neg.bed
	cat $? > $@

Fernandez_blood_cpg_pos.bed:	Supplemental_Table_4.txt
	join -t "	" -1 1 -2 4 <(sed '1,2d' $< | cut -f1 | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Fernandez_blood_cpg_neg.bed:	Supplemental_Table_5.txt
	join -t "	" -1 1 -2 4 <(sed '1,2d' $< | cut -f1 | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Fernandez_blood_cpg_all.bed:	Fernandez_blood_cpg_pos.bed Fernandez_blood_cpg_neg.bed
	cat $? > $@

