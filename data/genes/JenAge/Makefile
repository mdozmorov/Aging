all:	genes_JenAge.bed

genes_JenAge.bed:	gene.csv
	sed '1d' gene.csv | cut -d "," -f 24,25 | tr -d '\"' |awk -F "," '{if ($$2 == "Homo sapiens") print $$1}'| sort | uniq > genes_JenAge.txt && \
	python refgene.py genes_JenAge.txt | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@