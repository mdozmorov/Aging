#all:	Hannum_aging_markers_primary.bed Hannum_aging_markers_all.bed
#all:	Hannum_aging_markers_breast.bed Hannum_aging_markers_kidney.bed Hannum_aging_markers_lung.bed
all:	Hannum_genes_associated.bed Hannum_genes_predictive.bed Hannum_genes_predictive_pos.bed Hannum_genes_predictive_neg.bed

Hannum_aging_markers_primary.bed:	mmc2_Model_Primary_data.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	echo "chr2	30561970	30562093	ch.2.30415474F" >> $@ ; \
	echo "chr13	40666907	40667030	ch.13.39564907R" >> $@ ; \
	rm tmp.txt

Hannum_aging_markers_all.bed:	mmc2_Model_All_data.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	echo "chr13	40666907	40667030	ch.13.39564907R" >> $@ ; \
	rm tmp.txt

Hannum_aging_markers_breast.bed:	mmc3_Model_Breast.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Hannum_aging_markers_kidney.bed:	mmc3_Model_Kidney.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Hannum_aging_markers_lung.bed:	mmc3_Model_Lung.txt
	sed '1d' $< | cut -f1 | sort | uniq > tmp.txt && \
	join -t "	" -1 1 -2 4 tmp.txt ../../Illumina/Illumina450K.bed | awk 'BEGIN {OFS="\t"} {print $$2,$$3,$$4,$$1,$$5,$$6}' > $@ && \
	rm tmp.txt

Hannum_genes_associated.bed:	mmc4_Genes_S6.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Hannum_genes_predictive.bed:	mmc5_Genes_all_S7.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Hannum_genes_predictive_pos.bed:	mmc5_Genes_pos_S7.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@

Hannum_genes_predictive_neg.bed:	mmc5_Genes_neg_S7.txt
	python refgene.py $< | bedtools sort | mergeBed -s -c 4 -o distinct -i - > $@
