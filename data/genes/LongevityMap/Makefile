all:	genes_LongevityMap.bed

genes_LongevityMap.bed:	longevity.txt
	sed '1d' $< | cut -f6 | tr -d "\"" | tr "," "\n" | sort | uniq > genes_LongevityMap.txt && \
	python refgene.py genes_LongevityMap.txt | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@