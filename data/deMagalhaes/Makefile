#all:	deMagalhaes_genes_over.bed

all:	deMagalhaes_genes_under.bed

deMagalhaes_genes_over.bed:	deMagalhaes_genes_over.txt
	sed '1,15d' $< | cut -f3 | sort | uniq > genes_deMagalhaes_over.txt && \
	python refgene.py genes_deMagalhaes_over.txt | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

deMagalhaes_genes_under.bed:	deMagalhaes_genes_under.txt
	sed '1,10d' $< | cut -f3 | sort | uniq > genes_deMagalhaes_under.txt && \
	python refgene.py genes_deMagalhaes_under.txt | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@
