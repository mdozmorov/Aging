all:	Rakyan_aging_markers_all.bed

SHELL=/bin/bash

Rakyan_aging_markers_all.bed:	Supp_Table_3.txt
	join -t "	" -1 4 -2 1 <(cat ../../Illumina/Illumina27K.hg19.bed) <(cat $< | sed '1d' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@
