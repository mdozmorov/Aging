all:	Rakyan_cpg_all.bed

SHELL=/bin/bash

Rakyan_cpg_all.bed:	Supp_Table_3.txt
	join -t "	" -1 4 -2 1 <(sort -k4,4 ../../Illumina/Illumina450K.bed) <(cat $< | sed '1d' | cut -f1 | sort) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@
