all:	Hannum_cpg_primary.bed Hannum_cpg_all.bed
#all:	Hannum_cpg_breast.bed Hannum_cpg_kidney.bed Hannum_cpg_lung.bed
#all:	Hannum_genes_associated.bed Hannum_genes_predictive.bed Hannum_genes_predictive_pos.bed Hannum_genes_predictive_neg.bed

SHELL=/bin/bash

Hannum_cpg_primary.bed:	mmc2_Model_Primary_data.txt
	join -t "	" -1 1 -2 4 <(sed '1d' $< | cut -f1 | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Hannum_cpg_all.bed:	mmc2_Model_All_data.txt
	join -t "	" -1 1 -2 4 <(sed '1d' $< | cut -f1 | sort | uniq) <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@

Hannum_cpg_breast.bed:	mmc3_Model_Breast.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Hannum_cpg_kidney.bed:	mmc3_Model_Kidney.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Hannum_cpg_lung.bed:	mmc3_Model_Lung.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt <(sort -k4,4 ../../Illumina/Illumina450K.bed) | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Hannum_genes_associated.bed:	mmc4_Genes_S6.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Hannum_genes_predictive.bed:	mmc5_Genes_all_S7.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Hannum_genes_predictive_pos.bed:	mmc5_Genes_pos_S7.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Hannum_genes_predictive_neg.bed:	mmc5_Genes_neg_S7.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@
