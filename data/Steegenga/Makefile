all:	Steegenga_cpg_and_gene.bed Steegenga_cpg_not_gene.bed Steegenga_cpg_meta.bed
#all:	Steegenga_gene_and_cpg.bed
#all:	Steegenga_gene_not_cpg.bed

SHELL=/bin/bash

ILMN="/Users/mikhail/Documents/Work/GenomeRunner/Aging/Illumina/Illumina450K.bed"

Steegenga_cpg_and_gene.bed:	Steegenga_DMR_and_gene.txt
	join -t "	" -1 1 -2 4 <(sort -k1,1 $<) <(sort -k4,4 $(ILMN)) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Steegenga_cpg_not_gene.bed:	Steegenga_DMR_not_gene.txt
	join -t "	" -1 1 -2 4 <(sort -k1,1 $<) <(sort -k4,4 $(ILMN)) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Steegenga_cpg_meta.bed:	Steegenga_DMR_meta.txt
	join -t "	" -1 1 -2 4 <(sort -k1,1 $<) <(sort -k4,4 $(ILMN)) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Steegenga_gene_and_cpg.bed:	Steegenga_gene_and_DMR.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Steegenga_gene_not_cpg.bed:	Steegenga_gene_not_DMR.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@