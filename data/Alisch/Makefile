all:	Alisch_aging_markers_all.bed Alisch_aging_markers_pos.bed Alisch_aging_markers_neg.bed

SHELL=/bin/bash

Alisch_aging_markers_all.bed:	Alisch_et_al_Sup_Table2.txt
	join -t "	" -1 4 -2 1 <(cat ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Alisch_aging_markers_pos.bed:	Alisch_et_al_Sup_Table2.txt
	join -t "	" -1 4 -2 1 <(cat ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | awk '{if ($$13 > 0) print $$0}' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Alisch_aging_markers_neg.bed:	Alisch_et_al_Sup_Table2.txt
	join -t "	" -1 4 -2 1 <(cat ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | awk '{if ($$13 < 0) print $$0}' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@
