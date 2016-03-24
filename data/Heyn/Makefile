all:	Heyn_cpg_pos.bed Heyn_cpg_neg.bed Heyn_cpg_all.bed

SHELL=/bin/bash

Heyn_cpg_pos.bed:	sd02.txt
	join -t "	" -1 1 -2 4 <(sed '1,4d' $< | awk '$$4 == "Hypermethylated" {print $$5}' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Heyn_cpg_neg.bed:	sd02.txt
	join -t "	" -1 1 -2 4 <(sed '1,4d' $< | awk '$$4 == "Hypomethylated" {print $$5}' | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Heyn_cpg_all.bed:	Heyn_cpg_pos.bed Heyn_cpg_neg.bed
	cat $? | bedtools sort -i - | uniq > $@