all:	Alisch_cpg_all.bed Alisch_cpg_pos.bed Alisch_cpg_neg.bed

SHELL=/bin/bash

Alisch_cpg_all.bed:	Alisch_et_al_Sup_Table2.txt
	join -t "	" -1 4 -2 1 <(sort -k4,4 ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Alisch_cpg_pos.bed:	Alisch_et_al_Sup_Table2.txt
	join -t "	" -1 4 -2 1 <(sort -k4,4 ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | awk '{if ($$13 > 0) print $$0}' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Alisch_cpg_neg.bed:	Alisch_et_al_Sup_Table2.txt
	join -t "	" -1 4 -2 1 <(sort -k4,4 ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | awk '{if ($$13 < 0) print $$0}' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@
