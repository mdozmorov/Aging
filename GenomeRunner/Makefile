all:	cell_counts_all.txt cell_counts_all_10.txt

# Creates a list of all counts, in all experiments. The resulting file needs to be manually edited to remove artefacts.
ENCCELLS="/Users/mikhail/Documents/Work/GenomeRunner/genomerunner_database/ENCODE_cells.txt"

SHELL=/bin/bash

cell_counts_all.txt:
	for file in `find . -type f -name cell_counts.txt`; do cat $$file | cut -f1-3 >> $@; done && \
	sort -k2,2n -r -o $@ $@

cell_counts_all_10.txt:
	for file in `find . -type f -name cell_counts_10.txt`; do cat $$file | cut -f1-3 >> tmp.txt; done && \
	join -t "	" -a 1 <(sort tmp.txt) <(sort $(ENCCELLS)) > tmp1.txt && \
	sort -k2,2n -r tmp1.txt > $@ && \
	rm tmp*

# To do once, simply change the folder-specific Makefiles to count signal with higher threshold
change_p_val:
	for file in `find . -type f -name Makefile`; do sed 's/\ >\ 2/\ >\ 10/g' $$file | sed 's/cell\_counts/cell\_counts\_10/g' > `dirname $$file`"/Makefile10" ; done

rerun_p_val:
	for file in `find . -mindepth 1 -name Makefile10`; do cd `dirname $file`; make -f Makefile10; cd -; done

