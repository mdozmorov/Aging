all:	genes_GenAge.bed

genes_GenAge.bed:	genage_human.csv
	sed '1d' $< | cut -d "," -f3 | sort | uniq > genes_GenAge.txt && \
	python refgene.py genes_GenAge.txt | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@